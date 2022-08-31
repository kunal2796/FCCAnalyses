#python examples/FCCee/KG/preSel.py

from config.common_defaults import deffccdicts
import os

basedir=os.path.join(os.getenv('FCCDICTSDIR', deffccdicts), '') + "yaml/FCCee/spring2021/IDEA/"
outdir="outputs/FCCee/KG/"

#fcc_dir = os.getcwd().split("FCCAnalyses",1)[0]
#inputana=fcc_dir+"FCCAnalyses/examples/FCCee/KG/analysis.py"
#inputana="afs/cern.ch/work/k/kgautam/private/latest/FCCAnalyses/examples/FCCee/KG/analysis.py"

import multiprocessing
NUM_CPUS = int(multiprocessing.cpu_count()-2)

process_list=['p8_ee_Zuds_ecm91']
fraction=0.000099 #Zuds@91:100K(Zbb&Zcc)
#fraction=0.00099 #Zuds@91:1M(Zbb&Zcc)
#fraction=0.099  #Zqq@240
#fraction=0.016  #Zqq@365

import config.runDataFrame as rdf
myana=rdf.runDataFrame(basedir,process_list)
myana.run(ncpu=NUM_CPUS,fraction=fraction,outDir=outdir)
#myana.run(ncpu=NUM_CPUS,fraction=fraction,outDir=outdir,inputana=inputana)
