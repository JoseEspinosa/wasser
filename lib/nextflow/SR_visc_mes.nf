/*
#################################################################################
### Jose Espinosa-Carrasco. CB/CSN-CRG. Febraury 2016                         ###
#################################################################################
### Code : 02.03                                                              ### 
### SR pipeline                                                               ###
#################################################################################
*/

path_files = "$HOME/2015_viscMes/data/SR_analysis/"
params.in_trajectories_annot = "traj_annotated.csv"
params.dict_trans_found = "dict_traj_found.json"
params.lib_path="/git/wasser/"
traj_path = "${path_files}${params.in_trajectories_annot}"
dict_trans_path =  "${path_files}${params.dict_trans_found}"

viscMes_file = file(traj_path)
transitions_file = file(dict_trans_path)


R_lib="$HOME${params.lib_path}lib/R/"

//\'first\',\'second\'
sex_opt = ['m', 'w', 'm w']
//sex_opt = ['m']
//sex_opt = ['m', 'w']
//age_opt = Channel.from('30 35', '45')
// There are any men with an age of 70 or 75
age_opt = ['30 35', '40 45', '50 55', '60 65', '60 65 70 75 80', '30 35 40 45 50 55 60 65 70 75 80']
//age_opt = ['30 35']
//age_opt = ['30']

SR_mode =  ['counts', 'dosage']
//SR_mode =  ['counts']

process combine {
  input:
//  file bed_features.first()
//  file traj_annotation from viscMes_file.first()
//  file transitions_dict from transitions_file.first()
  file traj_annotation from viscMes_file
  file transitions_dict from transitions_file
  val sex from sex_opt
  each age from age_opt
  each mode from SR_mode
  
  //$HOME/2015_viscMes/lib/python/trans_bySubset_SR.py -i $HOME/2015_viscMes/data/SR_analysis/traj_annotated.csv -t $HOME/2015_viscMes/data/SR_analysis/dict_traj_found.json -s 'm' -a 30 -m 'dosage'
  output: 
  set file('*_tbl.csv'), sex, age, mode into sr_by_group
  set file('*_tbl.csv') into sr_by_group2write
  set file ('*.pdf') into heatmaps_SR
  script:
  println "Options for SR calculation are: $traj_annotation $transitions_dict $sex $age $mode"
  
//  each color from 'red','blue'
//  each size from 1,2
  """
  $HOME/git/wasser/lib/python/trans_bySubset_SR.py -i $traj_annotation -t $transitions_dict -s \'$sex\' -a \'$age\' -m $mode
  """
}

//Saving tables with transitions
sr_by_group2write.subscribe { 
    println "Received: " + it.name
    it.copyTo (it.name) 
}

//Saving heatmaps with transitions
heatmaps_SR.subscribe {
    println "Received: " + it.name
    it.copyTo (it.name)
}

//atc_df<- read.csv(file="/Users/jespinosa/Downloads/ATC.csv",header=TRUE, sep=",")
params.ATC_codes = "ATC.csv"

atc_path = "${path_files}${params.ATC_codes}"
atc_code_file = file(atc_path)

process top_20_hits {
    input:
    set file ('tbl_hits'), val (sex), val (age), val (mode) from sr_by_group
    set file ('ATC_ontology') from atc_code_file
    output:
    set file ('*.csv') into top_hits    

    script:
    joined_sex = sex.iterator().join('').replaceAll("\\s", "_")
    joined_age = age.iterator().join('').replaceAll("\\s", "_")
    joined_mode = mode.iterator().join('').replaceAll("\\s", "_")
    joined_options = joined_sex + "_" + joined_age + "_" + joined_mode

    println "Options for SR calculation are: $joined_options"
    println "Files input tables: $atc_code_file"
    
    """
    export R_LIBS="/software/R/packages"
    R --vanilla --slave --args ${tbl_hits} ${ATC_ontology} ${joined_options} < ${R_lib}top_hits_SR_names.R
    """        
}

//Saving tables with transitions
top_hits.subscribe {
    println "Received: " + it.name
    it.copyTo (it.name)
}


//export R_LIBS="/software/R/packages"
        
//R --vanilla --slave --args $tbl_hits $ATC_ontology  < ${R_lib}hits_annotation.R
   
/*
process heatmap_R {
  input: 
  file mat, sex, age, mode into sr_by_group
  
  output: 
  set file('*.png') into heat_map_R
  
}
*/
//export R_LIBS="/software/R/packages"


// For each subset calcuate SR both for the group selected (sex and age) 
/*process viscMes_subset {

 input:
 file 'viscMes_data' from viscMes_file
 file 'traj_data' from traj_file

 output:
 file '*.CULO' into subsets_viscMes


  """
  $HOME/2015_viscMes/lib/python/annotate_visc_tbl.py traj_data viscMes_data
  """
}

*/
// Old way of performing the script
//params.for_traj = "/Users/jespinosa/2015_viscMes/data/SR_analysis/for_traj.csv"
//params.rev_traj = "/Users/jespinosa/2015_viscMes/data/SR_analysis/rev_traj.csv"
//params.list_drugs = "/Users/jespinosa/2015_viscMes/data/SR_analysis/list_drugs.csv" 

// Combine trajectories with info from id to 
//process int2binFile {
//
// input:
// file 'traj_data' from traj_file
// file 'viscMes_data' from viscMes_file
// 
// output:
// file '*.csv' into subsets_viscMes

//
//  """
  //$HOME/2015_viscMes/lib/python/annotate_visc_tbl.py traj_data viscMes_data
//  """
//}


//para cada m w y edad hacer
// el archvio ya gurdado y hacer todos los filtros asi tengo ya los label



// I need flatten so it takes each item of the list one by one
// Here I get the cage from the file name and I put it on the channel
//bin_fileSingleCage = bin_fileSingle.flatten().map { binSingleCage -> 
//  def pattern = binSingleCage.name =~/^cage(\d+).*$/
//  println binSingleCage.name
//  println pattern [0][1]
//  def cage = pattern[0][1]
//  [ cage, binSingleCage ]
//}
