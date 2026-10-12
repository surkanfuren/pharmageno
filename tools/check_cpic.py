"""Inspect/update the same public CPIC store used by server.py (stdlib only)."""
import argparse
import json
import os
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from cpic_resources import ResourceStore

parser = argparse.ArgumentParser(description=__doc__)
actions = parser.add_mutually_exclusive_group()
actions.add_argument('--check', action='store_true', help='Force two complete live reads and stage any update')
actions.add_argument('--activate', metavar='SNAPSHOT_ID', help='Activate an explicitly reviewed staged update')
actions.add_argument('--passport', type=Path, help='Report impact on a saved passport without changing it')
args = parser.parse_args()
root = os.environ.get('CPIC_DATA_DIR', str(Path(__file__).resolve().parents[1] / '.cpic-data'))
store = ResourceStore(root)
if args.check:
    result = store.check(force=True)
elif args.activate:
    result = store.activate(args.activate)
elif args.passport:
    result = store.passport_impact(json.loads(args.passport.read_text()))
else:
    result = store.status()
print(json.dumps(result, indent=2))
if result.get('error') or result.get('freshness') in ('unverified', 'uninitialized'):
    sys.exit(1)
