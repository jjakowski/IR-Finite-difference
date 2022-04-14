#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <math.h>
using namespace std;

int main()
{
  double T1[1000];
  double tiny=1.e-16;
  double Pi=3.1415926535897932385;
  int i,j,n;
  char file1[50];
  FILE *inp1,*inp2,*out1;
  printf(" Which file ? \n");
  scanf("%s",file1);
  inp1 = fopen(file1,"rb");
  printf(" How many DP values to read ? \n");
  scanf("%d",&n);
  i=fread(T1,8,n,inp1);
//  for(i=0;i<n;i++){cout<<T1[i]<<" ";}cout<<"\n";
  for(i=0;i<n;i++){cout<< std::setprecision (15)  <<T1[i]<<" ";}cout<<"\n";
  return(0);
}
