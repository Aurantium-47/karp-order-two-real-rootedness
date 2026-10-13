"""Check file identity only; do not run or import mathematical verification code."""
from pathlib import Path
import hashlib
import json

root = Path(__file__).resolve().parent
manifest = json.loads((root / 'SHA256SUMS.json').read_text(encoding='utf-8'))
for item in manifest['files']:
    path = (root / item['path']).resolve()
    path.relative_to(root)
    data = path.read_bytes()
    if len(data) != item['bytes'] or hashlib.sha256(data).hexdigest() != item['sha256']:
        raise SystemExit('File mismatch: ' + item['path'])
print(json.dumps({'status': 'PASS', 'files_checked': len(manifest['files']),
                  'scope': 'File integrity only; no mathematical calculation.'}))
