//hello.cpp
#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <sstream>
#include <boost/algorithm/string.hpp>
using namespace std;

class PdbReader
{
private:
    string filename;
    int model_num;
public:
    vector<string> MODEL;
    PdbReader(string, int);
    int fileread();
    int movecenter();
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

int PdbReader::movecenter(){
    for(string line : MODEL){
        float x, y, z;
        std::vector<std::string> strs;
        boost::split(strs, line, boost::is_space(), boost::algorithm::token_compress_on);
        // strs is splitted_list of a line by spaces
    }
}

int main(){
    PdbReader a("2DX3.pdb", 4);
    a.fileread();
}