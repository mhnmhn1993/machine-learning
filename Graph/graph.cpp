#include "graph.h"

Graph::Graph(std::string name)
{
    std::ifstream in(filename);
    if (!in.is_open()){
        std::cout << "Error opening file"<<std::endl;
        return;
    }
    for(int i=1;!in.eof();i++){
        std::string line;
        getline(in,line);
        //analyze_adj(line,i);
    }
}

bool Graph::read_header(int& size, int& chain1 , int& consume1, int &cost1)
{
    std::ifstream in(filename);
    if(!in.is_open()){
        std::cout<<"Error opening file"<<std::endl;
        return false;
    }
    std::string line;
    getline(in,line);


    int count=0;
    int flag1=0;
    int flag2=0;
    int flag3=0;
    std::string node="";
    std::string chain="";
    std::string consume="";
    std::string cost="";
    int receiving=0;
    for(int i=0;!in.eof();i++){
        if(receiving==0 && line[i]==' '){
            continue;
        }
        if(flag1==0 && line[i]<='9' && line[i]>='0'){
            node+=line[i];
            receiving=1;
            continue;
        }
        if(flag1==0 && receiving==1 && line[i]==' '){
            node+='\0';
            receiving=0;
            size=atoi(node.c_str());
            flag1=1;
            count++;
            continue;
        }

        if(flag1==1 && flag2==0 && line[i]<='9' && line[i]>='0'){
            chain+=line[i];
            receiving=1;
            continue;
        }
        if(flag1==1 && flag2==0 && receiving==1 && line[i]==' '){
            chain+='\0';
            receiving=0;
            chain1=atoi(chain.c_str());
            flag2=1;
            count++;
            continue;
        }

        if(flag1==1 && flag2==1 && flag3==0 && line[i]<='9' && line[i]>='0'){
            consume+=line[i];
            receiving=1;
            continue;
        }
        /*if(flag1==1 && flag2==1 && flag3==0 && receiving==1 && line[i]==' '){
            consume+='\0';
            receiving=0;
            consume1=atoi(consume.c_str());
            flag3=1;
            count++;
            continue;
        }*/
    }
    consume+='\0';
    count++;
    consume1=atoi(consume.c_str());

    if(count!=3) return false;
    getline(in,cost);
    cost1=atoi(cost.c_str());
    return true;
}

int Graph::getNode()
{
    return node_size;
}

int Graph::getChain()
{
    return chain_size;
}

int Graph::getConsume()
{
    return consume_size;
}

int Graph::getCost()
{
    return cost;
}

/*bool Graph::analyze_adj(std::string line, int i)
{

}*/
