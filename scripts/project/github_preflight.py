"""Inspect Git-eligible files without staging, committing, pushing or scientific runs.

Reports contain paths, hashes and redacted finding categories, never matched secrets.
Every candidate is hashed; text is scanned, PNG structure is checked. This is a
bounded preflight, not a guarantee of privacy, licensing or absence of secrets.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import re
import struct
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
REPORT_REL = 'docs/audit/2026-10-27/github-sync'
TEXT_EXT = {'.r', '.py', '.sh', '.ps1', '.md', '.json', '.yml', '.yaml', '.toml',
            '.ini', '.cfg', '.tsv', '.csv', '.txt', '.svg', '.mol', '.sdf', '.smi', '.html'}
ALLOWED_DIRS = {'analysis', 'config', 'data', 'docs', 'environments', 'results', 'scripts', 'tests'}
ROOT_FILES = {'.gitignore', '.gitattributes', 'README.md', 'AGENTS.md', 'requirements.txt',
              '参考文献_方法依据.md', '单细胞基础补全_结果报告.md', '第一阶段_分析方案与验收.md',
              '第一阶段_结果报告.md', '工作流3_交接.md', '新对话接手_当前状态.md',
              '乳腺癌47靶点_查新检索记录.html'}
SECRET_PATTERNS = {
    'private_key': re.compile(r'-----BEGIN [A-Z ]*PRIVATE KEY-----'),
    'github_token': re.compile(r'\b(?:gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{30,})\b'),
    'openai_token': re.compile(r'\bsk-(?:proj-)?[A-Za-z0-9_-]{24,}\b'),
    'aws_access_key': re.compile(r'\b(?:AKIA|ASIA)[A-Z0-9]{16}\b'),
    'credential_literal': re.compile(r'''(?i)(?:api[_-]?key|access[_-]?token|auth[_-]?token|client[_-]?secret|secret[_-]?access[_-]?key|password)\s*[=:]\s*["']([^"'\n]{20,})["']'''),
    'credential_in_url': re.compile(r'https?://[^\s/:@]+:[^\s/@]+@'),
}
IDENTIFIERS = re.compile(r'(?i)^(patient_name|full_name|first_name|last_name|email|phone|telephone|date_of_birth|dob|address|medical_record_number|mrn)$')
PUBLIC_CASE_FIELDS = re.compile(r'(?i)(patient|donor|subject|case_id|clinical|treatment|subtype|ihc|pam50|age|sex|gender)')


def git(*args):
    result = subprocess.run(['git', *args], cwd=ROOT, capture_output=True, check=True)
    return result.stdout


def inspect_file(rel):
    path = ROOT / rel
    findings = []
    if path.is_symlink() or not path.resolve().is_relative_to(ROOT):
        return None, [dict(path=rel, line=0, kind='symlink_or_escape', severity='block')]
    payload = path.read_bytes()
    size = len(payload)
    first = Path(rel).parts[0]
    if first not in ALLOWED_DIRS and rel not in ROOT_FILES:
        findings.append(dict(path=rel, line=0, kind='outside_allowlist', severity='block'))
    if size > 10 * 1024**2:
        findings.append(dict(path=rel, line=0, kind='over_10MiB_requires_separate_review', severity='block'))
    suffix = path.suffix.lower()
    classification = 'text_or_small_result'
    if suffix == '.png':
        classification = 'png_figure'
        if not payload.startswith(b'\x89PNG\r\n\x1a\n'):
            findings.append(dict(path=rel, line=0, kind='invalid_png', severity='block'))
        else:
            offset = 8
            while offset + 12 <= len(payload):
                n = struct.unpack('>I', payload[offset:offset+4])[0]
                kind = payload[offset+4:offset+8]
                if kind in (b'tEXt', b'zTXt', b'iTXt', b'eXIf'):
                    findings.append(dict(path=rel, line=0, kind='png_metadata_'+kind.decode(), severity='review'))
                offset += n + 12
    elif suffix in TEXT_EXT or rel in ('.gitignore', '.gitattributes'):
        try:
            text = payload.decode('utf-8-sig')
        except UnicodeError:
            findings.append(dict(path=rel, line=0, kind='non_utf8_text', severity='block'))
            text = ''
        if '\0' in text:
            findings.append(dict(path=rel, line=0, kind='nul_in_text', severity='block'))
        for number, line in enumerate(text.splitlines(), 1):
            for kind, pattern in SECRET_PATTERNS.items():
                if pattern.search(line):
                    findings.append(dict(path=rel, line=number, kind=kind, severity='block'))
        if suffix in ('.csv', '.tsv') and text:
            header = next(csv.reader([text.splitlines()[0]], delimiter='\t' if suffix == '.tsv' else ','))
            if any(IDENTIFIERS.match(h.strip()) for h in header):
                findings.append(dict(path=rel, line=1, kind='direct_identifier_header', severity='block'))
            if any(PUBLIC_CASE_FIELDS.search(h) for h in header):
                classification = 'case_or_clinical_fields_review_source'
        if suffix in ('.md', '.html', '.svg', '.json') and re.search(r'(?i)(?:data:|file:///|src=["\'](?:https?:|[A-Z]:))', text):
            findings.append(dict(path=rel, line=0, kind='embedded_or_external_reference', severity='review'))
    else:
        findings.append(dict(path=rel, line=0, kind='unapproved_extension', severity='block'))
    return dict(path=rel, bytes=size, worktree_sha256=hashlib.sha256(payload).hexdigest(), classification=classification), findings


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check-reports', action='store_true', help='inspect existing generated reports without rewriting them')
    args = parser.parse_args()
    if Path(git('rev-parse', '--show-toplevel').decode().strip()).resolve() != ROOT:
        raise RuntimeError('Git root differs from project root')
    if args.check_reports:
        paths = sorted(p.relative_to(ROOT).as_posix() for p in (ROOT/REPORT_REL).glob('*') if p.is_file())
    else:
        paths = set(git('ls-files', '-z', '--cached', '--others', '--exclude-standard').decode('utf-8').strip('\0').split('\0'))
        paths = sorted(p for p in paths if p and not p.startswith(REPORT_REL + '/'))
    candidates, flags = [], []
    for rel in paths:
        row, findings = inspect_file(rel)
        if row:
            candidates.append(row)
        flags.extend(findings)
    samples = {'原始数据/example.tsv': True, 'outputs/example.txt': True,
               'data/external/example.tsv': True, 'outside.bam': True, 'outside.fastq.gz': True,
               'outside.h5ad': True, 'private.pem': True, '.env': True,
               'README.md': False, 'docs/TASK_LOGS/TEMPLATE.md': False,
               'results/10_programs/frozen/Programs_v1.json': False}
    for name, expected in samples.items():
        result = subprocess.run(['git', 'check-ignore', '--no-index', '-q', '--', name], cwd=ROOT)
        if result.returncode not in (0, 1) or (result.returncode == 0) != expected:
            raise RuntimeError('Ignore rule mismatch: ' + name)
    summary = dict(files=len(candidates), bytes=sum(r['bytes'] for r in candidates),
                   largest_bytes=max((r['bytes'] for r in candidates), default=0),
                   blocked=sum(f['severity'] == 'block' for f in flags),
                   manual_review_flags=sum(f['severity'] == 'review' for f in flags),
                   case_or_clinical_tables=sum(r['classification'] == 'case_or_clinical_fields_review_source' for r in candidates),
                   ignore_tests=len(samples), scientific_analysis_run=False, git_mutation=False,
                   limitation='Pattern scan is not a complete privacy, secret or licensing review; manually verify flags and clinical sources.')
    if args.check_reports:
        print(json.dumps(dict(mode='reports_read_only', **summary), ensure_ascii=False, indent=2))
        return 2 if summary['blocked'] else 0
    out = ROOT / REPORT_REL
    out.mkdir(parents=True, exist_ok=True)
    for name, rows, fields in [('UPLOAD_CANDIDATES.tsv', candidates, ['path', 'bytes', 'worktree_sha256', 'classification']),
                               ('REVIEW_FLAGS.tsv', flags, ['path', 'line', 'kind', 'severity'])]:
        with (out/name).open('w', encoding='utf-8', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t', lineterminator='\n')
            writer.writeheader()
            writer.writerows(rows)
    (out/'PREFLIGHT.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    local = out/'local'
    local.mkdir(exist_ok=True)
    (local/'candidate_paths.nul').write_bytes(b'\0'.join(r['path'].encode('utf-8') for r in candidates)+b'\0')
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 2 if summary['blocked'] else 0


if __name__ == '__main__':
    sys.exit(main())
