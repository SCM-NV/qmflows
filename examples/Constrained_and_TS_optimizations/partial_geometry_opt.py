from noodles import gather
from qmworks import adf, run, Settings, templates, molkit

methanol = molkit.from_smiles('CO')

# optimize only H atoms
s = Settings()
s.freeze = [0, 1]
geom_job1 = adf(templates.geometry.overlay(s), methanol, job_name='geom_job1').molecule

# optimize only H atoms
s = Settings()
s.selected_atoms = ['H']
geom_job2 = adf(templates.geometry.overlay(s), methanol, job_name='geom_job2').molecule

geom1, geom2 = run(gather(geom_job1, geom_job2), n_processes=1)

print(geom1)
print(geom2)