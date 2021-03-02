#include <iostream>
#include <fstream>
#include </usr/include/eigen3/Eigen/Dense>
#include <string>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

int mat2csv(MatrixXd mat, string filename)
{
      std::ofstream myfile;
      myfile.open (filename);
      for(int line = 0; line < mat.rows(); line++)
      {
	for(int column = 0; line < mat.cols(); column++)
	  {
	    myfile << to_string(line) + "," + to_string(mat(line,column)) + ",";
	  }
	myfile << "\n";
      }
      myfile.close();
      return 0;
}

int vec2csv(VectorXd vec, string filename)
{
      std::ofstream myfile;
      myfile.open (filename);
      for(int line = 0; line < vec.size(); line++)
      {
        myfile << to_string(line) + "," + to_string(vec(line)) + "\n";
      }
      myfile.close();
      return 0;
}
