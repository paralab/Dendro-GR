#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>

std::vector<std::string> rot_patterns;
std::vector<std::string> unique_rot_patterns;

#define SWAP(a, b) ((&(a) == &(b)) || \
                    (((a) -= (b)), ((b) += (a)), ((a) = (b) - (a))))


void rotateM(int index,int* current,int * rot_index,int dim)
{

  if(dim==2)
  {
    //index=rot_index[index];
    if(index==0)
    {
      rot_index[current[1]]=3;
      rot_index[current[3]]=1;
      SWAP(current[1],current[3]); // RIGHT Rotate and flip orientation


    }else if (index==1)
    {
      rot_index[current[0]]=2;
      rot_index[current[2]]=0;
      SWAP(current[0],current[2]); //LEFT Rotate and flip orientation:

    }

  }else if(dim==3)
  {
    //index=rot_index[index];
    if(index==0)
    {
      rot_index[current[1]]=7;
      rot_index[current[7]]=1;
      SWAP(current[1],current[7]);
      rot_index[current[4]]=2;
      rot_index[current[2]]=4;
      SWAP(current[2],current[4]);

    }else if(index==2)
    {
      rot_index[current[3]]=7;
      rot_index[current[7]]=3;
      SWAP(current[3],current[7]);
      rot_index[current[2]]=6;
      rot_index[current[6]]=2;
      SWAP(current[2],current[6]);

    }else if(index==1)
    {
      rot_index[current[3]]=5;
      rot_index[current[5]]=3;
      SWAP(current[3],current[5]);

      rot_index[current[3]]=7;
      rot_index[current[7]]=3;
      SWAP(current[3],current[7]);

      rot_index[current[2]]=6;
      rot_index[current[6]]=2;
      SWAP(current[2],current[6]);

      rot_index[current[0]]=2;
      rot_index[current[2]]=0;
      SWAP(current[0],current[2]);
    }else if(index==5)
    {
      rot_index[current[1]]=7;
      rot_index[current[7]]=1;
      SWAP(current[1],current[7]);

      rot_index[current[1]]=5;
      rot_index[current[5]]=1;
      SWAP(current[1],current[5]);

      rot_index[current[0]]=4;
      rot_index[current[4]]=0;
      SWAP(current[0],current[4]);

      rot_index[current[0]]=2;
      rot_index[current[2]]=0;
      SWAP(current[0],current[2]);

    }else if(index==6)
    {
      rot_index[current[1]]=5;
      rot_index[current[5]]=1;
      SWAP(current[1],current[5]);

      rot_index[current[0]]=4;
      rot_index[current[4]]=0;
      SWAP(current[0],current[4]);

    }else if(index==4)
    {

      rot_index[current[0]]=6;
      rot_index[current[6]]=0;
      SWAP(current[0],current[6]);

      rot_index[current[3]]=5;
      rot_index[current[5]]=3;
      SWAP(current[3],current[5]);
    }

  }

}


void rotate(int index,int* current,int * rot_index,int dim)
{

  if(dim==2)
  {
    //index=rot_index[index];
    if(index==0)
    {
      rot_index[current[1]]=3;
      rot_index[current[3]]=1;
      SWAP(current[1],current[3]); // RIGHT Rotate and flip orientation


    }else if (index==3)
    {
      rot_index[current[0]]=2;
      rot_index[current[2]]=0;
      SWAP(current[0],current[2]); //LEFT Rotate and flip orientation:

    }

  }else if(dim==3)
  {

    //index=rot_index[index];
    if(index==0)
    {
      rot_index[current[1]]=7;
      rot_index[current[7]]=1;
      SWAP(current[1],current[7]);
      rot_index[current[4]]=2;
      rot_index[current[2]]=4;
      SWAP(current[2],current[4]);

    }else if(index==1)
    {
      rot_index[current[3]]=7;
      rot_index[current[7]]=3;
      SWAP(current[3],current[7]);
      rot_index[current[2]]=6;
      rot_index[current[6]]=2;
      SWAP(current[2],current[6]);

    }else if(index==3)
    {
      rot_index[current[3]]=5;
      rot_index[current[5]]=3;
      SWAP(current[3],current[5]);

      rot_index[current[3]]=7;
      rot_index[current[7]]=3;
      SWAP(current[3],current[7]);

      rot_index[current[2]]=6;
      rot_index[current[6]]=2;
      SWAP(current[2],current[6]);

      rot_index[current[0]]=2;
      rot_index[current[2]]=0;
      SWAP(current[0],current[2]);
    }else if(index==4)
    {
      rot_index[current[1]]=7;
      rot_index[current[7]]=1;
      SWAP(current[1],current[7]);

      rot_index[current[1]]=5;
      rot_index[current[5]]=1;
      SWAP(current[1],current[5]);

      rot_index[current[0]]=4;
      rot_index[current[4]]=0;
      SWAP(current[0],current[4]);

      rot_index[current[0]]=2;
      rot_index[current[2]]=0;
      SWAP(current[0],current[2]);

    }else if(index==6)
    {
      rot_index[current[1]]=5;
      rot_index[current[5]]=1;
      SWAP(current[1],current[5]);

      rot_index[current[0]]=4;
      rot_index[current[4]]=0;
      SWAP(current[0],current[4]);

    }else if(index==7)
    {

      rot_index[current[0]]=6;
      rot_index[current[6]]=0;
      SWAP(current[0],current[6]);

      rot_index[current[3]]=5;
      rot_index[current[5]]=3;
      SWAP(current[3],current[5]);
    }

  }

}

void insert_unique(std::string str)
{
  for (int i=0;i<unique_rot_patterns.size();i++)
  {
    if(unique_rot_patterns[i].compare(str)==0){
      return;
    }
  }
  //std::cout<<"New Rotation:"<<str<<std::endl;
  rot_patterns.push_back(str);
  unique_rot_patterns.push_back(str);
}


int hilbert_rotation_pattern_2d()
{
  int default_rot []={0,2,3,1};
  //insert_unique("02315764");
  insert_unique("0231");
  int rot_count=0;
  int rot[4][4];
  int rot_current[4];
  int rot_index[4];
  while(rot_patterns.size()!=0)
  {
    std::string rotation=rot_patterns[0];
    rot_patterns.erase(rot_patterns.begin());
    for(int i=0;i<4;i++)
    {
      rot_current[i]=rotation.at(i)-'0';
    }

    for(int i=0;i<4;i++){
      rot[i][0]=rot_current[0];
      rot[i][1]=rot_current[1];
      rot[i][2]=rot_current[2];
      rot[i][3]=rot_current[3];

      rotate(i,rot[i],rot_index,2);
      char str[4]; //buffer size so never overflow
      sprintf( str, "%d%d%d%d", rot[i][0],rot[i][1],rot[i][2],rot[i][3]);
      //std::cout<<"new rotation: "<<str<<std::endl;
      insert_unique(std::string(str));

    }

  }


  return 0;
}



void generateTable_2d()
{

  //std::sort(unique_rot_patterns.begin(),unique_rot_patterns.end());
  //std::sort(unique_rot_patterns.begin(),unique_rot_patterns.end());
  char rot_index1[5];
  std::string rot_index_str;
  std::string rot_pat_index;
  for(int i=0;i<unique_rot_patterns.size();i++)
  {
    rot_index1[0]= unique_rot_patterns[i].find('0')+'0';
    rot_index1[1]= unique_rot_patterns[i].find('1')+'0';
    rot_index1[2]= unique_rot_patterns[i].find('2')+'0';
    rot_index1[3]= unique_rot_patterns[i].find('3')+'0';
    rot_index1[4]= '\0';
    rot_index_str=std::string(rot_index1);
    rot_pat_index=unique_rot_patterns[i]+rot_index_str;

    //std::cout<<"strcpy(rotations + "<<i*8<<", \""<<rot_pat_index<<"\");"<<std::endl;

    int rot[4][4];
    int rot_current[4];
    int rot_index[4];

    std::string rotation=unique_rot_patterns[i];
    //std::cout<<"Rotation:"<<rotation<<std::endl;
    //rot_patterns.erase(rot_patterns.begin());
    for(int w=0;w<4;w++)
    {
      rot_current[w]=rotation.at(w)-'0';
    }

    for(int w=0;w<4;w++){

      rot[w][0]=rot_current[0];
      rot[w][1]=rot_current[1];
      rot[w][2]=rot_current[2];
      rot[w][3]=rot_current[3];

      rotate(rot_index_str.at(w)-'0',rot[w],rot_index,2);
      char str[4]; //buffer size so never overflow
      sprintf( str, "%d%d%d%d", rot[w][0],rot[w][1],rot[w][2],rot[w][3]);
      //insert_unique(std::string(str));
      std::string rot_str=std::string(str);
      //std::cout<<"rotation"<<rot_str<<std::endl;
      int found=-1;
      for(int m=0;m<unique_rot_patterns.size();m++)
      {

        if(rot_str.compare(unique_rot_patterns[m])==0)
        {
          found=m;
          break;
        }
      }

      printf("HILBERT_TABLE[%d] = %d; \n" ,(4*i+w),found);

    }




  }


}

int hilbert_rotation_pattern()
{

  std::cout<<"Starting h_rotation tables"<<std::endl;

  int default_rot []={0,3,1,2,7,4,6,5};
  //insert_unique("02315764");
  insert_unique("01234567");
  int rot_count=0;
  int rot[8][8];
  int rot_current[8];
  int rot_index[8];
  while(rot_patterns.size()!=0)
  {
    std::string rotation=rot_patterns[0];
    rot_patterns.erase(rot_patterns.begin());
    for(int i=0;i<8;i++)
    {
      rot_current[i]=rotation.at(i)-'0';
    }

    for(int i=0;i<8;i++){
      rot[i][0]=rot_current[0];
      rot[i][1]=rot_current[1];
      rot[i][2]=rot_current[2];
      rot[i][3]=rot_current[3];
      rot[i][4]=rot_current[4];
      rot[i][5]=rot_current[5];
      rot[i][6]=rot_current[6];
      rot[i][7]=rot_current[7];

      rotate(i,rot[i],rot_index,3);
      char str[8]; //buffer size so never overflow
      sprintf( str, "%d%d%d%d%d%d%d%d", rot[i][0],rot[i][1],rot[i][2],rot[i][3],rot[i][4],rot[i][5],rot[i][6],rot[i][7]);
      //std::cout<<"Rotation pattern found: "<<str<<std::endl;
      insert_unique(std::string(str));

    }

  }
  return 0;
}



void generateTable()
{

  //std::sort(unique_rot_patterns.begin(),unique_rot_patterns.end());
  //std::sort(unique_rot_patterns.begin(),unique_rot_patterns.end());
  char rot_index1[9];
  std::string rot_index_str;
  std::string rot_pat_index;
  for(int i=0;i<unique_rot_patterns.size();i++)
  {
    rot_index1[0]= unique_rot_patterns[i].find('0')+'0';
    rot_index1[1]= unique_rot_patterns[i].find('1')+'0';
    rot_index1[2]= unique_rot_patterns[i].find('2')+'0';
    rot_index1[3]= unique_rot_patterns[i].find('3')+'0';
    rot_index1[4]= unique_rot_patterns[i].find('4')+'0';
    rot_index1[5]= unique_rot_patterns[i].find('5')+'0';
    rot_index1[6]= unique_rot_patterns[i].find('6')+'0';
    rot_index1[7]= unique_rot_patterns[i].find('7')+'0';
    rot_index1[8]='\0';
    rot_index_str=std::string(rot_index1);
    rot_pat_index=unique_rot_patterns[i]+rot_index_str;

    //std::cout<<"strcpy(rotations + "<<i*16<<", \""<<rot_pat_index<<"\");"<<std::endl;

    int rot[8][8];
    int rot_current[8];
    int rot_index[8];

    std::string rotation=unique_rot_patterns[i];
    //std::cout<<"Rotation:"<<rotation<<std::endl;
    //rot_patterns.erase(rot_patterns.begin());
    for(int w=0;w<8;w++)
    {
      rot_current[w]=rotation.at(w)-'0';
    }

    for(int w=0;w<8;w++){

      rot[w][0]=rot_current[0];
      rot[w][1]=rot_current[1];
      rot[w][2]=rot_current[2];
      rot[w][3]=rot_current[3];
      rot[w][4]=rot_current[4];
      rot[w][5]=rot_current[5];
      rot[w][6]=rot_current[6];
      rot[w][7]=rot_current[7];

      rotate(rot_index_str.at(w)-'0',rot[w],rot_index,3);
      char str[8]; //buffer size so never overflow
      sprintf( str, "%d%d%d%d%d%d%d%d", rot[w][0],rot[w][1],rot[w][2],rot[w][3],rot[w][4],rot[w][5],rot[w][6],rot[w][7]);
      //insert_unique(std::string(str));
      std::string rot_str=std::string(str);
      //std::cout<<"rotation"<<rot_str<<std::endl;
      int found=-1;
      for(int m=0;m<unique_rot_patterns.size();m++)
      {

        if(rot_str.compare(unique_rot_patterns[m])==0)
        {
          found=m;
          break;
        }
      }

      printf("HILBERT_TABLE[%d] = %d; \n" ,(8*i+w),found);


    }




  }


}


int main(int argc, char **argv) {
  //std::cout << "Hello, world!" << std::endl;
  //hilbert_rotation_pattern();
  //generateTable();

  hilbert_rotation_pattern_2d();
  generateTable_2d();

  return 0;
}
