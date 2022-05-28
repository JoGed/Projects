#ifndef OPTFUNCTIONS_H
#define OPTFUNCTIONS_H

void ReadAndSendParams();
void CheckInput();
void InitializeOR();
void InitializeDR();
void SetupORDRCommunication();
void CleanUpDRInit();
void ZeroDensInOR(int r);
void InitializeSimParams();
void SaveRegions();
void AdjSource(int r);
void UpdateGradients(int r);
void StoreEfieldDR(int r);
void StoreEfieldORAndComputeObj(int r);
void CheckObjChange();
REAL SendObjToOptRanks();
void SetupForwardSim();
void SetupAdjointSim();
void UpdateDesign();
void MeanDFT(int r);
void AllocTopOptFields();
void ReleaseTopOptFields();
void FinalizeOpt();

#endif