#Creation of job_scripts for SLURM queue handler

import os
import shutil
import pickle
import re
from pydantic import BaseModel, parse_obj_as
from openfe.setup.methods.openmm.equil_rbfe_methods import (
    SystemSettings, TopologySettings, AlchemicalSettings,
    OpenMMEngineSettings, SamplerSettings, BarostatSettings,
    IntegratorSettings, SimulationSettings
)
from openff.units import unit

class SLURM():

    class SLURMSettings(BaseModel):
        job_name = "FEP_model_job"
        account = "def-schapira"
        ntasks = "4"
        mem_per_cpu = "2G"
        mail_usr = ""
        mail_type = ""
        time = "0-10:00:00"

    class SLURMSettingsM(BaseModel):
        job_name = "FEP_master_job"
        account = "def-schapira"
        ntasks = "32"
        mem_per_cpu = "2G"
        mail_usr = ""
        mail_type = ""
        time = "0-10:00:00"
    
    def write_job_script(SLURMSettings, directory, filename="job_script.sh"):
        write_file = open(filename,'w')
        write_file.write("#!/bin/bash\n")
        write_file.write("#SBATCH --job-name="+SLURMSettings.job_name+"\n")
        write_file.write("#SBATCH --account="+SLURMSettings.account+"\n")
        write_file.write("#SBATCH --ntasks="+SLURMSettings.ntasks+"\n")
        write_file.write("#SBATCH --mem-per-cpu="+SLURMSettings.mem_per_cpu+"\n")
        if (SLURMSettings.mail_usr!="") & (SLURMSettings.mail_type !=""):
            write_file.write("#SBATCH --mail-user="+SLURMSettings.mail_usr+"\n")
            write_file.write("#SBATCH --mail-type="+SLURMSettings.mail_type+"\n")
        write_file.write("#SBATCH --time="+SLURMSettings.time+"\n")
        write_file.write("\n")
        write_file.write("cd "+directory+"\n")
        write_file.write("conda activate openfe\n")
        write_file.write("python/run.py complex_transform.pkl solvent_transform.pkl\n")
        write_file.write("while true\n")
        write_file.write("do\n")
        write_file.write("    if [[ ! -f 'solvent_transform_final.pkl' ]]\n")
        write_file.write("    then\n")
        write_file.write("        sleep 120\n")
        write_file.write("        continue\n")
        write_file.write("    else\n")
        write_file.write("        break\n")
        write_file.write("    fi\n")
        write_file.write("done\n")
        write_file.write("conda deactivate\n")
        write_file2 = open("run.py",'w')
        write_file2.write("import os\n")
        write_file2.write("import sys\n")
        write_file2.write("import copy\n")
        write_file2.write("import pickle\n")
        write_file2.write("import from openfe.setup import ChemicalSystem\n")
        write_file2.write("from openfe.setup import SolventComponent, ProteinComponent\n")
        write_file2.write("from openfe.setup import Network\n")
        write_file2.write("from openfe.setup import SmallMoleculeComponent\n")
        write_file2.write("from openfe.setup.ligand_network_planning import generate_radial_network\n")
        write_file2.write("from openfe.setup.lomap_mapper import LomapAtomMapper\n")
        write_file2.write("from openfe.setup.methods.openmm.equil_rbfe_methods import (SystemSettings, TopologySettings, AlchemicalSettings, OpenMMEngineSettings, SamplerSettings, BarostatSettings, IntegratorSettings, SimulationSettings)\n")
        write_file2.write("from openfe.setup.methods.openmm.equil_rbfe_methods import RelativeLigandTransformSettings\n")
        write_file2.write("from openfe.setup.methods.openmm.equil_rbfe_methods import RelativeLigandTransform\n")
        write_file2.write("from openfe.setup.slurmscripts import SLURM\n")
        write_file2.write("\n")
        write_file2.write("open_object(sys.arg[1])\n")
        write_file2.write("open_object(sys.arg[2])\n")
        write_file2.write("complex_transform.run(verbose=True)\n")
        write_file2.write("solvent_transform.run(verbose=True)\n")
        write_file2.write("SLURM.save_object(complex_transform, 'complex_transform_final.pkl')\n")
        write_file2.write("SLURM.save_object(solvent_transform, 'solvent_transform_final.pkl')\n")

    def write_master_script(SLURMSettingsM, ligand_file, filename="job_master_script.sh", nof_jobs=0):
        if (nof_jobs==0):
            print("Number of jobs is not specified")
        write_file = open(filename,'w')
        write_file.write("#!/bin/bash\n")
        write_file.write("#SBATCH --job-name="+SLURMSettingsM.job_name+"\n")
        write_file.write("#SBATCH --account="+SLURMSettingsM.account+"\n")
        write_file.write("#SBATCH --ntasks="+SLURMSettingsM.ntasks+"\n")
        write_file.write("#SBATCH --mem-per-cpu="+SLURMSettingsM.mem_per_cpu+"\n")
        if (SLURMSettingsM.mail_usr!="") & (SLURMSettingsM.mail_type !=""):
            write_file.write("#SBATCH --mail-user="+SLURMSettingsM.mail_usr+"\n")
            write_file.write("#SBATCH --mail-type="+SLURMSettingsM.mail_type+"\n")
        write_file.write("#SBATCH --time="+SLURMSettingsM.time+"\n")
        write_file.write("\n")
        write_file.write("for i in {1.."+str(nof_jobs)+"..1}\n")
        write_file.write("do\n")
        write_file.write("    while true\n")
        write_file.write("    do\n")
        write_file.write("        a=$(sq -h -r | grep 'FEP_model_job' | wc -l)\n")
        write_file.write("        if (( $a < 8 ))\n")
        write_file.write("        then\n")
        write_file.write("            job_name = './molcule_'+$i+'/job_script.sh'\n")
        write_file.write("            srun $job_name &\n")
        write_file.write("            break\n")
        write_file.write("        else\n")
        write_file.write("            sleep 120\n")
        write_file.write("            continue\n")
        write_file.write("        fi\n")
        write_file.write("    done\n")
        write_file.write("done\n")
        write_file.write("wait\n")
        write_file.write("\n")
        write_file.write("python analysis.py "+str(nof_jobs)+"\n")
        write_file2 = open("analysis.py",'w')
        write_file2.write("import os\n")
        write_file2.write("import sys\n")
        write_file2.write("import copy\n")
        write_file2.write("import pickle\n")
        write_file2.write("from rdkit.Chem import PandasTools\n")
        write_file2.write("import from openfe.setup import ChemicalSystem\n")
        write_file2.write("from openfe.setup import SolventComponent, ProteinComponent\n")
        write_file2.write("from openfe.setup import Network\n")
        write_file2.write("from openfe.setup import SmallMoleculeComponent\n")
        write_file2.write("from openfe.setup.ligand_network_planning import generate_radial_network\n")
        write_file2.write("from openfe.setup.lomap_mapper import LomapAtomMapper\n")
        write_file2.write("from openfe.setup.methods.openmm.equil_rbfe_methods import (SystemSettings, TopologySettings, AlchemicalSettings, OpenMMEngineSettings, SamplerSettings, BarostatSettings, IntegratorSettings, SimulationSettings)\n")
        write_file2.write("from openfe.setup.methods.openmm.equil_rbfe_methods import RelativeLigandTransformSettings\n")
        write_file2.write("from openfe.setup.methods.openmm.equil_rbfe_methods import RelativeLigandTransform\n")
        write_file2.write("from openfe.setup.slurmscripts import SLURM\n")
        write_file2.write("df = PandasTools.LoadSDF('"+ligand_file+"', embedProps=True, molColName=None, smilesName='smiles')\n")
        write_file2.write("dG_complex = []\n")
        write_file2.write("error_complex = []\n")
        write_file2.write("dG_solvent = []\n")
        write_file2.write("error_solvent = []\n")
        write_file2.write("for i in range(1,"+str(nof_jobs)+")\n")
        write_file2.write("    open_object('./molecule'+str(i)+'/complex_transform_final.pkl'\n")
        write_file2.write("    open_object('./molecule'+str(i)+'/solvent_transform_final.pkl'\n")
        write_file2.write("    complex_results = complex_transform.get_results()\n")
        write_file2.write("    solvent_results = solvent_transform.get_results()\n")
        write_file2.write("    dG_complex.append(complex_results.dG())\n")
        write_file2.write("    error_complex.append(complex_results.dG_error())\n")
        write_file2.write("    dG_solvent.append(solvent_results.dG())\n")
        write_file2.write("    error_solvent.append(solvent_results.dG_error())\n")
        write_file2.write("df['dG_complex'] = dG_complex\n")
        write_file2.write("df['error_complex']=error_complex\n")
        write_file2.write("df['dG_solvent']=dG_solvent\n")
        write_file2.write("df['error_solvent']=error_solvent\n")

    def save_object(obj, filename):
        with open(filename, 'wb') as outp:  # Overwrites any existing file.
            pickle.dump(obj, outp, pickle.HIGHEST_PROTOCOL)

    def open_object(filename):
        with open(filename, 'wb') as inp:
            pickle.load(inp)
    def read_params(config_file: str) -> dict:
        with open(config_file) as f:
            lines = [l.strip() for l in f.readlines() if l.strip() != '']
        params = {k.strip(): v.strip() for k, v in [l.split(';') for l in lines]}
        return params
    def updated_params(params):
        systemsettings = SystemSettings.parse_raw(pickle.dumps(eval(params['systemsettings'])), content_type='application/pickle', allow_pickle=True)
        topologysettings = TopologySettings.parse_raw(pickle.dumps(eval(params['topologysettings'])), content_type='application/pickle', allow_pickle=True)
        alchemicalsettings = AlchemicalSettings.parse_raw(pickle.dumps(eval(params['alchemicalsettings'])), content_type='application/pickle', allow_pickle=True)
        samplersettings = SamplerSettings.parse_raw(pickle.dumps(eval(params['samplersettings'])), content_type='application/pickle', allow_pickle=True)
        barostatsettings = BarostatSettings.parse_raw(pickle.dumps(eval(params['barostatsettings'])), content_type='application/pickle', allow_pickle=True)
        integratorsettings = IntegratorSettings.parse_raw(pickle.dumps(eval(params['integratorsettings'])), content_type='application/pickle', allow_pickle=True)
        simulationsettings = SimulationSettings.parse_raw(pickle.dumps(eval(params['simulationsettings'])), content_type='application/pickle', allow_pickle=True)
        return systemsettings,topologysettings,alchemicalsettings,samplersettings,barostatsettings,integratorsettings,simulationsettings
        

