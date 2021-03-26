TString arrayOfPeriods[] = {"LHC16e2plus"};
int numberOfPeriods = sizeof(arrayOfPeriods) / sizeof(arrayOfPeriods[0]);

void CopyTree(const char *ext = ".root")
{

   for (int iPeriod = 0; iPeriod < numberOfPeriods; iPeriod++)
   {
      TString dirname;
      dirname.Form("%s/", arrayOfPeriods[iPeriod].Data());
      //Loop over files in subdirectory
      TSystemDirectory dir(dirname, dirname);
      TList *files = dir.GetListOfFiles();
      if (files)
      {
         TSystemFile *file;
         TString fname;
         TIter next(files);
         while ((file = (TSystemFile *)next()))
         {
            fname = file->GetName();
            if (!file->IsDirectory() && fname.EndsWith(ext))
            {
               cout << "Extracting tree from " << fname.Data() << endl;

               TString filename;
               filename.Form("%s/%s", dirname.Data(), fname.Data());

               TFile *oldfile = TFile::Open(filename);
               TObjArray *obj = ((TObjArray *)oldfile->Get("EventTrees_CAny"));
               TTree *oldtree = ((TTree *)obj->FindObject("eventsTree"));

               // Activate all branches
               oldtree->SetBranchStatus("*", 1);

               // Create a new file + a clone of old tree in new file
               TString newFileName = filename;
               newFileName.ReplaceAll("AnalysisResults", "AnalysisResultsNew");
               TFile *newfile = TFile::Open(newFileName, "recreate");
               auto newtree = oldtree->CloneTree();

               newtree->Print();
               newfile->Write();
               delete newfile;
               delete oldfile;
               //End of new tree
            }
         }
      }
      //end of loop over files in directory
      //////////////////////////////////////////////////
   }
}
