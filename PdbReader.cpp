//hello.cpp
#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <sstream>
#include <boost/algorithm/string.hpp>
#include <eigen3/Eigen/Dense>
using namespace std;
using namespace Eigen;
class PdbReader
{
private:
    string filename;
    int model_num;
public:
    vector<string> MODEL;
    PdbReader(string, int);
    int fileread();
    MatrixXd get_matrix(int, int);
    int file_changer(string additional_name);
};

PdbReader::PdbReader(string s, int num = -1){
    filename = s;
    model_num = num;
}
int PdbReader::fileread(){
    cout << filename << endl;
    std::ifstream ifs(filename);
    if (!ifs) {
        cerr << "Can't open file" << endl;
        return -1;
    }
    string line;
    vector< vector<string> > all_list;
    vector<std::string> each_list;
    if (model_num == -1){
        while (std::getline(ifs, line)){
            if (line.substr(0, 5) == "MODEL"){
                each_list.push_back(line);
            }
        }
        MODEL = each_list;
    }

    else{
        while (std::getline(ifs, line)){
            if (line.substr(0, 5) == "MODEL"){
                each_list.clear();
            }
            else if(line.substr(0, 4) == "ATOM"){
                each_list.push_back(line);
            }
            else if(line.substr(0, 6) == "ENDMDL"){
                all_list.push_back(each_list);
            }
        }
        MODEL = all_list[model_num - 1];
    }
    return 0;
}

MatrixXd PdbReader::get_matrix(int x_start_col, int element_col){
    int count(-1);
    int N(MODEL.size());
    int weight_sum;
    MatrixXd Q(N, 3);
    MatrixXd M(N, 3);
    for(string line : MODEL){
        count += 1;
        double x, y, z;
        std::vector<std::string> strs;
        // strs is splitted_list of a line by spaces
        boost::split(strs, line, boost::is_space(), boost::algorithm::token_compress_on);
        x = std::stod(strs[x_start_col - 1]);
        y = std::stod(strs[x_start_col]);
        z = std::stod(strs[x_start_col + 1]);
        Vector3d q(x, y, z);
        Q.row(count) = q;
        // calculate moment
        string element(strs[element_col -1]);
        int weight;
        if (element == "H"){
            weight = 1;
        }
        else if (element == "C"){
            weight = 12;
        }
        else if (element == "N"){
            weight = 14;
        }
        else if (element == "O"){
            weight = 16;
        }
        else if (element == "S"){
            weight = 32;
        }
        else{
            cout << "[WARNING] Write the weight of " + element << endl;
        }
        weight_sum += weight;
        RowVector3d moment;
        moment = weight * q;
        M.row(count) = moment;
    }
    RowVector3d c;
    MatrixXd Q_from_center(N, 3);
    c << M.col(0).sum()/weight_sum, M.col(1).sum()/weight_sum, 
    M.col(2).sum()/weight_sum;
    Q_from_center = Q.rowwise() - c;
    return Q_from_center;
}

int main(){
    MatrixXd Q1_from_center;
    PdbReader Reader1("2DX3.pdb", 4);
    Reader1.fileread();
    Q1_from_center = Reader1.get_matrix(7, 12);
    cout << Eigen::JacobiSVD(Q1_from_center) << endl;
}