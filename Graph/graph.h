#ifndef GRAPH_H
#define GRAPH_H
#include <string>
#include <fstream>
#include <iostream>
#include <istream>
class Graph{
private:
    std::string filename;
    int AdjMatrix[1000][1000];
    int node_size;
    int adj_flag=0;
    int chain_size;
    int consume_size;
    int cost;
public:
    Graph(std::string name);
    //~Graph();
    //void AdjMatrix();
    //bool analyze_adj(,std::string line ,int i);
    bool read_header(int& size, int& chain1 , int& consume1, int& cost1);
    int getNode();
    int getChain();
    int getConsume();
    int getCost();
};
















#endif // GRAPH_H
