import copy
import json
import random
import sys
import tempfile
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from core_balance import allocate_counts,fair,p25,LEAF
from reallocate_core import event_allocate,load_inputs

CFG=dict(balance=dict(foil='equal',delta=1.5,hole=1.5),qa_min=10,fit_max=200000,seed=667)


def fixture():
    rows=[]; data={}; meta={}
    for fi,f in enumerate([-8.,-3.,0.,3.,8.]):
        for setting in range(2):
            name=f'rg{fi}{setting}'
            meta[name]=dict(foils=[f],edges=[-10.,-8.,-5.])
            events=[]
            for nd,(lo,hi) in enumerate([(-10.,-8.),(-8.,-5.)]):
                for h,n in enumerate([1,5+fi,40+setting*10,200*(fi+1)]):
                    for i in range(n+3):
                        events.append(dict(entry=len(events)+2**54,foil=0,ndel=nd,zfoil=f,delta_low=lo,delta_high=hi,
                                           xscol=h,yscol=4,sample=2 if i==n else (0 if i>n else (1 if i%2 else 3)),quality=2 if i<=n else (i%2),core_score=0.9))
                    rows.append(dict(rungroup=name,zfoil=f,delta_low=lo,delta_high=hi,xscol=h,yscol=4,
                                     available=n,total=n+3,holdout=1,foil=0,ndel=nd))
            data[name]={k:np.array([r[k] for r in events]) for k in events[0]}
    return rows,data,dict(cfg=CFG,meta=meta,blocked={(2,3),(5,5)})


class BalanceTests(unittest.TestCase):
    def test_waterfill_and_stable_remainder(self):
        self.assertEqual(fair([100,300,600],700,['a','b','c'],667),[100,300,300])
        self.assertEqual(fair([100,300,600],450,['a','b','c'],667),[100,175,175])
        a=fair([3,3,3],2,['a','b','c'],667)
        self.assertEqual(a,fair([3,3,3],2,['c','b','a'],667)[::-1])
        self.assertEqual(sum(fair([10**12]*10,10**12,list(range(10)),667)),10**12)
    def test_p25(self):
        self.assertEqual(float(p25([0,1,100,200,400])),75.25)
        self.assertEqual(float(p25([])),0)
    def test_fixed_caps_and_exact_ids_reordering(self):
        rows,data,inp=fixture()
        r=allocate_counts(rows,CFG,[-8.,-3.,0.,3.,8.],[(-10.,-8.),(-8.,-5.)])
        before=copy.deepcopy(data)
        ids,_=event_allocate(data,inp,r)
        rows.reverse()
        data2={name:{k:v[::-1].copy() for k,v in a.items()} for name,a in reversed(list(before.items()))}
        r2=allocate_counts(rows,CFG,[-8.,-3.,0.,3.,8.],[(-8.,-5.),(-10.,-8.)])
        ids2,_=event_allocate(data2,inp,r2)
        self.assertEqual(ids,ids2);self.assertEqual(r['leaves'],r2['leaves'])
        self.assertEqual(len({r['quota'] for r in r['foils']}),1)
        for name,a in data.items():
            np.testing.assert_array_equal(a['quality'],before[name]['quality'])
            np.testing.assert_array_equal(a['sample']==2,before[name]['sample']==2)
            self.assertTrue(np.all(a['quality'][a['sample']==1]==2))
        for l in r['leaves']:
            self.assertTrue(0<=l['quota']<=l['capacity']<=l['available'])
        for d in r['deltas']:
            ls=[l for l in r['leaves'] if l['zfoil']==d['zfoil'] and l['delta_low']==d['delta_low']]
            self.assertEqual(sum(l['quota'] for l in ls),d['quota'])
            self.assertEqual(d['q_delta'],float(p25([sum(l['available'] for l in r['leaves'] if l['zfoil']==d['zfoil'] and l['delta_low']==lo) for lo in [-10.,-8.]])))
    def test_qa_reporting_only_and_zero_budget(self):
        rows,_,_=fixture()
        args=(rows,CFG,[-8.,-3.,0.,3.,8.],[(-10.,-8.),(-8.,-5.)])
        r=allocate_counts(*args)
        other=allocate_counts(rows,{**CFG,'qa_min':10000},args[2],args[3])
        self.assertEqual([x['quota'] for x in r['leaves']],[x['quota'] for x in other['leaves']])
        tiny=allocate_counts(rows,{**CFG,'balance':dict(foil='equal',hole=.001,delta=1.5)},args[2],args[3])
        self.assertEqual(tiny['training_total'],0);self.assertEqual(len(tiny['zero_capacity_foils']),5)
        small=allocate_counts(rows,{**CFG,'fit_max':7},args[2],args[3])
        self.assertEqual(small['training_total'],5)
        self.assertTrue(any(x['qa_low'] and x['quota']==0 for x in small['leaves']))
        empty=allocate_counts(rows,CFG,args[2]+[10.],args[3])
        self.assertEqual(empty['training_total'],0);self.assertIn(10.,empty['zero_capacity_foils'])
    def test_invalid(self):
        rows,_,_=fixture()
        for bad in (0,-1,float('nan'),float('inf'),True):
            with self.assertRaises(ValueError):allocate_counts(rows,{**CFG,'balance':dict(foil='equal',hole=bad,delta=1.5)},[-8.,-3.,0.,3.,8.],[(-10.,-8.),(-8.,-5.)])
        with self.assertRaises(ValueError):allocate_counts(rows+rows[:1],CFG,[-8.,-3.,0.,3.,8.],[(-10.,-8.),(-8.,-5.)])
    def test_blocked_holdout_and_counts(self):
        rows,data,inp=fixture()
        # Turn an entire observed label into a confirmed exclusion, including core holdout.
        inp['blocked'].add((0,4))
        for r in rows:
            if r['xscol']==0:r['available']=0
        result=allocate_counts(rows,CFG,[-8.,-3.,0.,3.,8.],[(-10.,-8.),(-8.,-5.)])
        ids,audit=event_allocate(data,inp,result)
        self.assertFalse(any(r['xscol']==0 for r in ids));self.assertTrue(audit)
        self.assertTrue(any(r['protected'] for r in audit))
        for a in data.values():self.assertFalse(np.any((a['sample']==1)&(a['balance_excluded']!=0)))
    def test_campaign_preview(self):
        campaign=Path(__file__).resolve().parents[1]/'HMS_6p667GeV'
        inp=load_inputs(campaign,'min10')
        result=allocate_counts(inp['rows'],inp['cfg'],inp['required'],zip(inp['edges'],inp['edges'][1:]))
        self.assertEqual(result['training_total'],68535)
        self.assertEqual(result['foil_budget'],13707)
        self.assertEqual(sum(r['available'] for r in result['leaves']),310513)
        self.assertEqual(len(result['deltas']),25)
        self.assertTrue(all(r['delta_loss']==0 for r in result['deltas']))

if __name__=='__main__':unittest.main()
