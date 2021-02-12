"""
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import argparse
try:
  import cPickle as pickle
except:
  import pickle
import json
import os
import sys

from numpy import array, random
import numpy as np
from scipy import stats

from . import aggregate, crystal, rotator, generator


if sys.version_info[0] >= 3:
    xrange = range


def generate_aggregate(monomer_generator,N=5,align=True):

    align_rot = rotator.PartialAligningRotator(exp_sig_deg=40)
    uniform_rot = rotator.UniformRotator()

    agg = [monomer_generator() for i in xrange(N)]

    while len(agg) > 1:
        r = array([((a.extent[0][1]-a.extent[0][0])+(a.extent[1][1]-a.extent[1][0]))/4.0 for a in agg])
        m_r = np.sqrt(array([a.X.shape[0] for a in agg])/r)
        r_mat = (np.tile(r,(len(agg),1)).T+r)**2
        mr_mat = abs(np.tile(m_r,(len(agg),1)).T - m_r)
        p_mat = r_mat * mr_mat
        p_mat /= p_mat.max()
        collision = False
        while not collision:
            
            i = random.randint(len(agg))
            j = random.randint(len(agg))
            rnd = random.rand()
            if rnd < p_mat[i][j]:
                print(i, j)
                agg_top = agg[i] if (m_r[i] > m_r[j]) else agg[j]
                agg_btm = agg[i] if (m_r[i] <= m_r[j]) else agg[j]
                agg_btm.rotate(uniform_rot)
                collision = agg_top.add_particle(particle=agg_btm.X,required=True,pen_depth=80e-6)
                if collision:
                    if align:
                        agg_top.align()
                        agg_top.rotate(align_rot)
                    else:
                        agg_top.rotate(uniform_rot)
                    agg.pop(i if (m_r[i] <= m_r[j]) else j)            

    if align:
        agg[0].align()
        agg[0].rotate(align_rot)
    agg[0].rotate(rotator.HorizontalRotator())

    return agg[0]


def gen_monomer(psd="monodisperse", size=1.0, min_size=1e-3, max_size=10,
    mono_type="dendrite", grid_res=0.02e-3, rimed=False):
        
    def make_cry(D):
        if mono_type=="dendrite":
           current_dir = os.path.dirname(os.path.realpath(__file__))
           with open(current_dir+"/dendrite_grid.dat", 'rb') as f:
               kwargs = {"encoding": "latin1"} if sys.version_info[0] >= 3 else {}
               grid = pickle.load(f, **kwargs)
           cry = crystal.Dendrite(D, hex_grid=grid)
        elif mono_type=="plate":
            cry = crystal.Plate(D)            
        elif mono_type=="needle":
            cry = crystal.Needle(D)
        elif mono_type=="rosette":
            cry = crystal.Rosette(D)
        elif mono_type=="bullet":
            cry = crystal.Bullet(D)
        elif mono_type=="spheroid":
            cry = crystal.Spheroid(D,0.6)
        return cry
                
    rot = rotator.UniformRotator()            
    
    def gen():
        if psd=="monodisperse":
            D = size
        elif psd=="exponential":
            psd_f = stats.expon(scale=size)
            D=max_size+1
            while (D<min_size) or (D>max_size):
                D = psd_f.rvs()
        
        cry = make_cry(D)
        
        gen = generator.MonodisperseGenerator(cry, rot, grid_res)
        if rimed:
            agg = aggregate.RimedAggregate(gen)
        else:
            agg = aggregate.Aggregate(gen)
        return agg
    
    return gen


def visualize_crystal(mono_type):
    gen = gen_monomer(mono_type=mono_type, size=2e-3, grid_res=40e-6)
    cry = gen()
    cry.align()
    cry.visualize(bgcolor=(1,1,1))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--psd', required=True)
    parser.add_argument('--mono_type', required=True)
    parser.add_argument('--mono_size', type=float, required=True)
    parser.add_argument('--mono_min_size', type=float, default=None)
    parser.add_argument('--mono_max_size', type=float, default=None)
    parser.add_argument('--num_monos', type=int, required=True)
    parser.add_argument('--output', type=argparse.FileType('w'), required=True)
    parser.add_argument('--grid_res', type=float, required=True)
    args = parser.parse_args()

    mono_generator = gen_monomer(psd=args.psd, size=args.mono_size, 
        min_size=args.mono_min_size, max_size=args.mono_max_size,
        mono_type=args.mono_type, grid_res=args.grid_res)
        
    agg = generate_aggregate(mono_generator,N=args.num_monos,align=True)

    meta = {"psd": args.psd, "mono_type": args.mono_type, 
        "mono_size": args.mono_size, "mono_min_size": args.mono_min_size,
        "mono_max_size": args.mono_max_size, "num_monos": args.num_monos,
        "grid_res": args.grid_res, "file_name": args.output.name,
        "extent": agg.extent}
    np.savetxt(args.output, agg.grid(), fmt="%d")
    json.dump(meta, file(args.output.name+".meta", 'w'))
