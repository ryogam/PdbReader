//hello.cpp
#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <sstream>
using namespace std;

class PdbReader
{
    string filename;
public:
    PdbReader(string s);
    int fileread();
    int atomout();
};

PdbReader::PdbReader(string s) : filename(s){}

int PdbReader::fileread(){
    std::ifstream ifs(filename);
    if (!ifs) {
        cerr << "Can't open file" << endl;
        return -1;
    }
    string line;
    vector< vector<string> > all_list;
    vector<std::string> each_list;
    while (std::getline(ifs, line)){
        if (line.substr(0, 5) == "MODEL"){
            all_list.push_back(each_list);
            each_list.clear();
            string model_name;
            model_name = line;
            std::cout << model_name << std::endl;
            each_list.push_back(model_name);
        }
        else if (line.substr(0, 4) == "ATOM"){
            istringstream iss(line);
            string s;
            while (iss >> s) {
                each_list.push_back(s);
            }
        }
    }
    return 0;
}


int main()
{
    PdbReader a("2DX3.pdb");
    a.fileread();
}