"""Read-only LFS object availability; authentication stays in process memory."""
from pathlib import Path
import subprocess,json,urllib.request,urllib.parse,urllib.error,argparse
from datetime import datetime,timezone

HERE=Path(__file__).resolve().parent
parser=argparse.ArgumentParser()
parser.add_argument('--repo',type=Path,default=HERE.parent if HERE.name=='scripts' else HERE/'github-pulizia-grafene')
parser.add_argument('--report',type=Path,default=None if HERE.name=='scripts' else HERE/'lfs-publication-readback.json')
args=parser.parse_args();repo=args.repo.resolve()
listed=json.loads(subprocess.run(['git','lfs','ls-files','--json'],cwd=repo,capture_output=True,text=True,check=True).stdout)
entries=listed['files'] if isinstance(listed,dict) else listed
objects=[dict(oid=x.get('oid') or x.get('sha256'),size=x['size']) for x in entries]
auth=subprocess.run(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=10','git@github.com',
                     'git-lfs-authenticate','mapice/pulizia-grafene.git','download'],
                    capture_output=True,text=True,timeout=20)
assert auth.returncode==0,'LFS download authentication failed'
access=json.loads(auth.stdout);url=access['href'].rstrip('/')+'/objects/batch'
host=urllib.parse.urlsplit(url).hostname;assert host=='github.com' or host.endswith('.github.com')
headers=dict(access.get('header',{}));headers.update({'Content-Type':'application/vnd.git-lfs+json',
                                                   'Accept':'application/vnd.git-lfs+json'})
request=urllib.request.Request(url,data=json.dumps(dict(operation='download',transfers=['basic'],objects=objects)).encode(),headers=headers)
with urllib.request.urlopen(request,timeout=30) as response:result=json.loads(response.read())
status=[]
for obj in result['objects']:
    action=obj.get('actions',{}).get('download');code=None;remote_length=None
    if action:
        address=action['href'];parsed=urllib.parse.urlsplit(address)
        assert parsed.scheme=='https'
        assert any(parsed.hostname==d or parsed.hostname.endswith('.'+d)
                   for d in ['github.com','githubusercontent.com','amazonaws.com'])
        h=dict(action.get('header',{}));h['Range']='bytes=0-0'
        try:
            with urllib.request.urlopen(urllib.request.Request(address,headers=h),timeout=20) as response:
                code=response.status
                content_range=response.headers.get('Content-Range')
                remote_length=int(content_range.rsplit('/',1)[-1]) if content_range else int(response.headers['Content-Length'])
                assert len(response.read(1))==1
        except urllib.error.HTTPError as error:code=error.code
    status.append(dict(oid=obj['oid'],bytes=obj['size'],download_offered=bool(action),
                       batch_error_code=obj.get('error',{}).get('code'),
                       range_response_status=code,remote_length_bytes=remote_length,
                       object_exists_with_expected_size=code in [200,206] and remote_length==obj['size']))
report=dict(time_UTC=datetime.now(timezone.utc).isoformat(),objects=status,
            all_objects_readback_available=all(r['object_exists_with_expected_size'] for r in status))
if args.report is not None:args.report.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps([dict(oid=r['oid'][:12],bytes=r['bytes'],available=r['object_exists_with_expected_size'],
                      status=r['range_response_status'] or r['batch_error_code']) for r in status]))
