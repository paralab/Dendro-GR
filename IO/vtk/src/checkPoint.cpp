//
// Created by milinda on 7/12/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief Contans utility functions to save/load mesh (octree) and variable list defined on the otree.
*/
//

#include "checkPoint.h"

namespace io
{

    namespace checkpoint
    {

      int writeOctToFile(const char * fName,const ot::TreeNode* pNodes,const unsigned int num)
      {
          FILE* outfile = fopen(fName,"w");
          if(outfile==NULL) {std::cout<<fName<<" file open failed "<<std::endl; return 1;}
          fwrite(&num,sizeof(unsigned int),1,outfile); // write out the number of elements.

          if(num>0)
            fwrite(pNodes,sizeof(ot::TreeNode),num,outfile);

          fclose(outfile);
          return 0;
      }


      int readOctFromFile(const char * fName,std::vector<ot::TreeNode> & pNodes)
      {
          FILE* inpfile = fopen(fName,"r");
          if(inpfile==NULL) {std::cout<<fName<<" file open failed "<<std::endl; return 1;}
          unsigned int num=0;
          fread(&num,sizeof(unsigned int ),1,inpfile);

          if(num>0)
          {
              pNodes.resize(num);
              fread(&(*(pNodes.begin())),(sizeof(ot::TreeNode)),num,inpfile);
          }

          fclose(inpfile);
          return 0;
      }


    }// end of namespace checkpoint

} // end of namespace io