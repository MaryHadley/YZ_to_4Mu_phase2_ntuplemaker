# YZ_to_4Mu_phase2_ntuplemaker

cd testBjorn  

#To run Z+Y->4 mu looper, do:   
root -l -b  
.L flexible_loop_big4MuVtxCutRemoved.C++    
run("<fileName.root>")   

#To run Z->4 mu looper, do:   
root -l -b   
.L flexible_loop_Zto4Mu.C++  
run("<fileName.root>")  

#If running loopers with January cuts omitted, use the clean_flexible versions! 
#Remember to set the triggerYear variable appropriately when running the flexible loopers!


