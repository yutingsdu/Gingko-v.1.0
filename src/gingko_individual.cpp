#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<algorithm>
#include<vector>
#include<map>
#include "gingko_individual.h"
using namespace std;
typedef vector<int> path_t;
ofstream outgtf;
int SampleSize;
double Coverage = 0;
double SEED = 0;
struct info
{
    int number;
    double seed;
    int normal_edge_number;
    int partial_edge_number;
    int mj_number;
    int Nmj_number;
    int seed_sample_number;
};
map<string, vector<double> > id_cov_map;
map<string, info > id_info_map;

void load_info(char*file)
{
    ifstream in(file);
    string s;
    istringstream istr;
    while(getline(in,s))
    {
        istr.str(s);
	string temp,id;
	double seed;
	int normal_edge,partial_edge,mj,Nmj,seed_sample_number;
	istr>>temp>>id>>seed>>normal_edge>>partial_edge>>temp>>temp>>temp>>mj>>Nmj>>seed_sample_number;
	istr.clear();
	int N = atoi(id.substr(id.length() - 3,1).c_str());
	info info_={N,seed,normal_edge,partial_edge,mj,Nmj,seed_sample_number};
	id_info_map[id] = info_;
    }
    //cout<<"id_info_map.size(): "<<id_info_map.size()<<endl;
}
void process(string tranid,vector<string> oneTrans)
{
    double r1,r2,r3,r4,r5;
    r1=0.0;r2=0.5;r3=1.5;r4=5.0;r5=10.0;
    /*
    if(SampleSize>30){
        r1=0.1;r2=1.5;r3=5.0;r4=10.0;r5=20.0;
    }
    else {
        r1=2.5;r2=1.0;r3=5.0;r4=10.0;r5=20.0;
    }*/
    if(id_info_map.find(tranid) != id_info_map.end() )//&& id_cov_map.find(tran_id) != id_cov_map.end())
    {
	    double cov = 0;
	    for(size_t i=0;i<id_cov_map[tranid].size();i++)
		    cov += id_cov_map[tranid][i];

	    bool flag = false;

	    info ti = id_info_map[tranid];
//	    if(SEED>=1 && ti.seed_sample_number < SampleSize) return;
//	    if(ti.Nmj_number >=1) return;
//	    if( ti.normal_edge_number<=1)return;
	    
	    if(id_info_map[tranid].number == 1)
	        if(cov >= r1*SEED) flag = true;

	    if(id_info_map[tranid].number == 2)
	        if(cov >= r2*SEED) flag = true;

	    if(id_info_map[tranid].number == 3)
	        if(cov > r3*SEED) flag = true;

	    if(id_info_map[tranid].number ==4 )
	        if(cov > r4*SEED) flag = true;

	    if(id_info_map[tranid].number >4 )
		    if(cov > r5*SEED) flag = true;
		    

	    if(flag)  
	    {
		    for(size_t i=0;i<oneTrans.size();i++)
		    outgtf<<oneTrans[i]<<endl;
	    }
     }
    return;
}
void get_final_results(char*file)
{

    ifstream in(file);
    istringstream istr;
    string s;

    string chr, strand, lable;
    string tranid;
    string temp, Cov_s;
    double Cov;
    int exon_l,exon_r;
    vector<int> vecExon;
    vector<string> oneTrans;
    getline(in,s);
    istr.str(s);
    while(istr>>temp)
	    if( temp == "transcript_id")
		    istr>>tranid;
    istr.clear();
    oneTrans.push_back(s);

    while(getline(in,s))
    {
	istr.str(s);
 	string current_id;
	string curr_chr, curr_strand;
	int curr_el, curr_er;

	istr>>curr_chr>>temp>>lable>>curr_el>>curr_er>>temp>>curr_strand;


	while(istr>>temp)
	{
	    if( temp == "transcript_id")
		istr>>current_id;
	}
	if(current_id ==  tranid)
	{
	  oneTrans.push_back(s);
	}
	else
	{
	  process(tranid,oneTrans);
	  tranid = current_id;
	  oneTrans.clear();
	  oneTrans.push_back(s);
	}

	istr.clear();
    }
    process(tranid,oneTrans);
    return;

}
//./exe a.gtf b.gtf .. N.gtf input.info input.gtf output.gtf DEED //only remove repeat ones
//
int main(int argc,char* argv[])
{
    //load_transref(argv[1],intron_trans_map,Chr_Junc_map);
    //cout<<argc<<endl;
    SampleSize = argc - 5;
    string S = argv[argc-1];
    SEED = atof(S.c_str());
    //cout<<"filter: "<<SEED<<endl;
    load_info(argv[argc-4]);

    outgtf.open(argv[argc-2]);

    for(int i=1;i<=argc-5;i++){
	//cout<<"load: "<<argv[i]<<" "<<i<<" sample..."<<endl;
        load_transref(argv[i],id_cov_map);
    
    }

    get_final_results(argv[argc-3]);
    outgtf.close();
    return 0;
}
