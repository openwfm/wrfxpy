from __future__ import absolute_import
from __future__ import print_function
from forecast import make_job_file, JobState, process_arguments, load_sys_cfg
import json
import logging
import sys,glob
import os.path as osp

simulations_path = osp.abspath('simulations')    # copy here job_id/wfc
catalog_path = osp.join(simulations_path,'catalog.json')
try:
    catalog = json.load(open(catalog_path,'r'))
except:
    print(('Cannot open catalog at %s, creating new.' % catalog_path))
    catalog={}

if __name__ == '__main__':

    if len(sys.argv) < 2:
        print('Usage: ./recover_catalog.sh 1.json 2.json ....')
        print('x.json are inputs to forecast.sh as from wrfxctrl, starting with /, no spaces')
        print('Example: ./recover_catalog.sh ~/Projects/wrfxctrl/jobs/*.json')
        print('Important: must be run from the wrfxpy directory. Before using:')
        print('In wrfxweb/fdds/simulations: tar cvfz ~/c.tgz <simulations to recove>/*.json catalog.json')
        print('Transfer the file c.tgz and  untar here')
        print('On exit, ./simulations/catalog.json will be updated')
        sys.exit(1)

    sims = glob.glob(osp.join(simulations_path,'wfc-*'))
    state = {}
    desc = {}
    for s in sims:
        job_id = osp.basename(s)
        if job_id in catalog:
            state[job_id]='Already in catalog'
            desc[job_id] = catalog[job_id]['description'] 
        else:
            state[job_id]='No catalog entry'
            desc[job_id] = 'Unknown' 
     
    for js_path in sys.argv[1:]:
        # print js_path
        js = json.load(open(js_path,'r'))
        description = js['postproc']['description']
        #print json.dumps(js, indent=4, separators=(',', ': '))
        jsb=osp.basename(js_path)
        m=glob.glob(osp.join(simulations_path,'*','wfc-'+jsb))
        print(('%s simulations found for %s file %s' % 
            (len(m), description, js_path)))
        for mf in m: 
            #print mf
            manifest=osp.basename(mf)
            job_id=osp.basename(osp.split(mf)[0])
            print('%s found in %s' % (manifest, job_id))
            manifest_js=json.load(open(mf,'r'))
            from_utc = '999-99-99_99:99:99'
            to_utc = '0000-00-00_00:00:00'
            for domain in manifest_js:
                for time in manifest_js[domain]:
                    from_utc = min(from_utc,time)
                    to_utc = max(to_utc,time)
            new={'manifest_path':osp.join(job_id,manifest),
                            'description':js['postproc']['description'],
                            'from_utc':from_utc,
                            'to_utc':to_utc}
            if job_id in catalog:
                if catalog[job_id]==new:
                    print('Catalog entry already exists, no change.')
                    state[job_id] = 'No change'
                else:
                    print(('Replacing catalog entry %s' % job_id))
                    print(('Old value:\n%s' % 
                        json.dumps(catalog[job_id], indent=4, separators=(',', ': '))))
                    print(('New value:\n%s' % 
                        json.dumps(new, indent=4, separators=(',', ': '))))
                    state[job_id] = 'Replaced'
            else:
                    print(('Creating catalog entry %s' % job_id))
                    print(('New value:\n%s' % 
                        json.dumps(new, indent=4, separators=(',', ': '))))
                    state[job_id] = 'Recovered'
            catalog[job_id]=new
            desc[job_id] = catalog[job_id]['description'] 
    print('Writing catalog at %s' % catalog_path)
    json.dump(catalog, open(catalog_path,'w'), indent=4, separators=(',', ': '))
    for job_id in sorted(state):
        print('%s: %s: %s' % (job_id, desc[job_id], state[job_id]))
