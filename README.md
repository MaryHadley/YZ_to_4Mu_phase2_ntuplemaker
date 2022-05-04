# YZ_to_4Mu_phase2_ntuplemaker

cd testBjorn  

#To run Z+Y->4 mu looper, do:   
root -l -b  
.L loop_big4MuVtxCutRemoved_12April2022.C++   
run( <fileName.root> )   

#To run Z->4 mu looper, do:   
root -l -b   
.L loop_Zto4Mu.C++   
run(<fileName.root>)  

