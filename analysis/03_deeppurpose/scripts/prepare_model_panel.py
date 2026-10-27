"""Download/inspect official linked archives; do not infer on the study compound."""
import concurrent.futures
import hashlib
import json
import pathlib
import pickle
import io
import time
import urllib.request
import zipfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
MODELS = {
    "CNN_CNN_BindingDB": 4159715, "Daylight_AAC_BindingDB": 4159667,
    "Morgan_AAC_BindingDB": 4159694, "Morgan_CNN_BindingDB": 4159708,
    "MPNN_CNN_BindingDB": 4204178,
    "CNN_CNN_BindingDB_IC50": 4203606, "Daylight_AAC_BindingDB_IC50": 4203604,
    "Morgan_AAC_BindingDB_IC50": 4203602, "Morgan_CNN_BindingDB_IC50": 4203603,
    "MPNN_CNN_BindingDB_IC50": 4203605,
}

class PlainConfigUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        raise pickle.UnpicklingError(f"Global class disallowed in configuration: {module}.{name}")

def prepare(item):
    name, file_id = item
    archive_dir = ROOT / "model_archives"
    archive_dir.mkdir(parents=True, exist_ok=True)
    path = archive_dir / (name + ".zip")
    url = f"https://dataverse.harvard.edu/api/access/datafile/{file_id}"
    record = {"model": name, "file_id_in_author_code": file_id, "source_url": url,
              "dataset_family": "BindingDB", "endpoint": "IC50" if name.endswith("_IC50") else "Kd",
              "architecture_selected_for_availability_validation_only": True}
    try:
        if not path.exists():
            if name == "Morgan_CNN_BindingDB":
                previous = ROOT.parent / "DeepPurpose_准备_2026-10-23" / "model_cache" / "model_morgan_cnn_bindingdb.zip"
                if previous.exists():
                    import shutil
                    shutil.copyfile(previous, path)
                    record["reused_previous_download"] = True
            if not path.exists():
                partial = path.with_suffix(".zip.part")
                for attempt in range(3):
                    try:
                        request = urllib.request.Request(url, headers={"User-Agent": "Workflow3-DeepPurpose-preparation/1.0"})
                        with urllib.request.urlopen(request, timeout=45) as response, partial.open("wb") as out:
                            length = response.headers.get("Content-Length")
                            total = 0
                            while chunk := response.read(1024 * 1024):
                                total += len(chunk)
                                if total > 300_000_000:
                                    raise ValueError("Archive exceeds 300 MB download limit")
                                out.write(chunk)
                            if length and total != int(length):
                                raise ValueError("Incomplete download")
                        with zipfile.ZipFile(partial) as z:
                            if z.testzip() is not None:
                                raise ValueError("ZIP CRC failure")
                        partial.replace(path)
                        break
                    except Exception:
                        if attempt == 2:
                            raise
                        time.sleep(2)
        with zipfile.ZipFile(path) as archive:
            if archive.testzip() is not None:
                raise ValueError("ZIP CRC failure")
            entries = archive.infolist()
            configs = [x for x in entries if pathlib.PurePosixPath(x.filename).name == "config.pkl" and not x.filename.startswith("__MACOSX/")]
            weights = [x for x in entries if pathlib.PurePosixPath(x.filename).name == "model.pt" and not x.filename.startswith("__MACOSX/")]
            if len(configs) != 1 or len(weights) != 1:
                raise ValueError("Expected exactly one configuration and one model.pt")
            config = PlainConfigUnpickler(io.BytesIO(archive.read(configs[0]))).load()
            destination = ROOT / "models" / name
            destination.mkdir(parents=True, exist_ok=True)
            for member, filename in [(configs[0], "config.pkl"), (weights[0], "model.pt")]:
                (destination / filename).write_bytes(archive.read(member))
            (destination / "config_review.json").write_text(json.dumps(config, indent=2), encoding="utf-8")
        record.update(status="archive_and_plain_configuration_verified_not_model_inference", bytes=path.stat().st_size,
                      archive_sha256=hashlib.sha256(path.read_bytes()).hexdigest(), config=config,
                      model_pt_sha256=hashlib.sha256((destination / "model.pt").read_bytes()).hexdigest())
    except Exception as exc:
        record.update(status="failed", error=str(exc))
    print(json.dumps({k: v for k, v in record.items() if k != "config"}, ensure_ascii=True), flush=True)
    return record

def main():
    (ROOT / "metadata").mkdir(parents=True, exist_ok=True)
    records = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        for record in pool.map(prepare, MODELS.items()):
            records.append(record)
            (ROOT / "metadata" / "Model_panel_archive_audit.json").write_text(json.dumps(records, ensure_ascii=False, indent=2), encoding="utf-8")
    if any(x["status"] == "failed" for x in records):
        raise SystemExit(1)

if __name__ == "__main__":
    main()
