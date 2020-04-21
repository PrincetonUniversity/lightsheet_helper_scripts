# Helper scripts to assist with Wang lab lightsheet & tracing projects & Princeton BRAIN CoGS Histology Core.

Written by Z. Dhanerawala and A. Hoag (zmd@princeton.edu, ahoag@princeton.edu).

### Relies on [BrainPipe](https://github.com/PrincetonUniversity/BrainPipe) and [ClearMapCluster](https://github.com/PrincetonUniversity/ClearMapCluster) local installations.

## Module functions:
* `projects` --> project-specific scripts written for Princeton Neuro researchers
* `clearmap` --> helper scripts and post-hoc processing for images analyzed using [ClearMapCluster](https://github.com/PrincetonUniversity/ClearMapCluster)
* `datajoint` --> scripts to injest and curate lightsheet data in DataJoint databases; related to A. Hoag lightserv project
* `lightsheet` --> related to processing and analyzing raw or downsampled lightsheet images 
* `neuroglancer` --> scripts to injest raw and downsampled lightsheet images into the Neuroglancer interface; related to A. Hoag lightserv project
* `ontology_analysis` --> scripts to assist in traversing an atlas hierarchy of structures; typically used for the Allen Mouse Brain Atlas
* `registration` --> helper scripts to do registration independent of the [BrainPipe](https://github.com/PrincetonUniversity/BrainPipe) pipeline
* `registration/annotation_volume_analysis` --> scripts utilizing atlas annotation volumes for processing or to do volume analysis

### Project specific module functions:
In `projects`:
* `braincogs` --> general analysis scripts written for BRAIN CoGS members/labs
* `esteban_engel` --> analysis scripts related to Engel lab LAP promoter project [Maturana et al. (2020)](https://www.sciencedirect.com/science/article/pii/S2329050120300644)
* `jv` --> cfos and general analysis scripts for J. Verpeut
* `malaz` --> analysis scripts to quantify clearing efficacy based on a summer undergraduate research project
* `mkislin` --> probe detection analysis scripts for M. Kislin in-vivo recording work
* `rat_atlas` --> scripts related to the Princeton Rat Atlas project in BRAIN CoGS
* `scooter` --> analysis scripts for Falkner lab's c-Fos projects
* `seagravesk` --> cfos and general analysis scripts for K. Seagraves
* `tracing` --> analysis for Wang lab tracing projects; includes the HSV and PRV scripts for Pisano et al. (2020)

