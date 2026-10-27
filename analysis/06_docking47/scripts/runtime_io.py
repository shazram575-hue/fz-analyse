"""Retry readers for a status/metadata file that may be inside a short write window."""
import json,time
from pathlib import Path
def read_json(path):
    path=Path(path)
    for i in range(20):
        try:return json.loads(path.read_text(encoding='utf8'))
        except (json.JSONDecodeError,PermissionError):
            if i==19:raise
            time.sleep(.1)
