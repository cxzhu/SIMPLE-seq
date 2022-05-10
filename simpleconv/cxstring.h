//
// Created by Chenxu Zhu on 6/29/18.
//

#ifndef REACHTOOLS_CXSTRING_H
#define REACHTOOLS_CXSTRING_H

#include <algorithm>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <deque>
#include <string>

using namespace std;

class cxstring {
private:

public:
    static vector<string> split(const string &s, const string &seperator);
    static string chomp(string str);
    static bool is_bam_header(const std::string& buffer);
    static int str2int(string str);
    static string int2str(int i);


};


class samline{
public:
    string readname, chr, cigar, chrnext, seq, qual, other;
    int flag, pos, mapq, posnext, plen;
        
    void init(string str){
        str = cxstring::chomp(str);
        vector<string> tmp = cxstring::split(str, "\t");
        readname = tmp[0];
        flag = cxstring::str2int(tmp[1]);
        chr = tmp[2];
        pos = cxstring::str2int(tmp[3]);
        mapq = cxstring::str2int(tmp[4]);
        cigar = tmp[5];
        chrnext = tmp[6];
        posnext = cxstring::str2int(tmp[7]);
        plen = cxstring::str2int(tmp[8]);
        seq = tmp[9];
        qual = tmp[10];
        if(tmp.size() > 11){
        	other = tmp[11];
        }
        else{
        	other = "";
        }
        return;
    }
    void empty(){
        readname = "";
        flag = 0;
        chr = "";
        pos = 0;
        mapq = 0;
        cigar = "";
        chrnext = "";
        posnext = 0;
        plen = 0;
        seq = "";
        qual = "";
        return;
    }
    void write(FILE * outfile){
    	string out = readname + "\t" + cxstring::int2str(flag) + "\t" + chr + "\t" + cxstring::int2str(pos) + "\t" +  cxstring::int2str(mapq) + "\t" + cigar + "\t" + chrnext +  "\t" + cxstring::int2str(posnext) + "\t" + cxstring::int2str(plen) + "\t" + seq + "\t" + qual + "\t" + other;
    	fputs((out+"\n").c_str(), outfile);
    }
};

class fqline{
public:
	string readname,seq,mark,qual;

	void read_part_record(FILE * infile, string rn){
		readname = rn;
		char buffer[2000];
		stringstream tmp;
		fgets(buffer, sizeof(buffer), infile); // seq
		tmp << buffer;
		tmp >> seq;
		tmp.str("");
		fgets(buffer, sizeof(buffer), infile); // mark
		tmp << buffer;
		tmp >> mark;
		tmp.str("");
		fgets(buffer, sizeof(buffer), infile); // qual
		tmp << buffer;
		tmp >> qual;
		return;
	}
	void read_full_record(FILE * infile){
		char buffer[2000];
		stringstream tmp;
		fgets(buffer, sizeof(buffer), infile); // readname
		tmp << buffer;
		tmp >> readname;
		tmp.str("");
		fgets(buffer, sizeof(buffer), infile); // seq
		tmp << buffer;
		tmp >> seq;
		tmp.str("");
		fgets(buffer, sizeof(buffer), infile); // mark
		tmp << buffer;
		tmp >> mark;
		tmp.str("");
		fgets(buffer, sizeof(buffer), infile); // qual
		tmp << buffer;
		tmp >> qual;
		return;
	}
	void write_record(FILE * outfile){
		fputs((readname + "\n").c_str(), outfile);
		fputs((seq + "\n").c_str(), outfile);
		fputs((mark + "\n").c_str(), outfile);
		fputs((qual + "\n").c_str(), outfile);
	}
};

class read2_2r {
private:
	int align_score(string str1, string str2){
		int score = 0;
    for(int i = 0; i < str1.length(); ++i){
      if(str2[i] == 'N')continue;
      str2[i] == str1[i] ? (score += 2) : (score -= 1);

    }
    return score;
	}
public:
	int bc1, bc2, bc4;
	string bc, rawline, sbc1, sbc2, sbc4, bsbc4, type, umi;
	int dock;
	bool valid;

	void init(string line){
		rawline = line;
		dock = -1;
		valid = false;
		return;
	}

	void trim(){

		//".   NNNNNNNNNNNNNNNNNGTGGCCGATGTTTCGGTGCGAACTCAGACCNNNNNNNATCCACGTGCTTGAGGCATTCGAGNNNNN";
		// 1 determine additional Ns

		if(rawline.length() < 93)return;
		dock = 0;
		int t = 0;
		int cur_s = 0;
		string bait = "GTGGCCGATGTTTCG";
		for(int i = -2; i < 7; ++i){
			string qu = rawline.substr(17+i, 15);
			int score = align_score(qu, bait);
			if(score < cur_s)continue;
			cur_s = score;
			t = i;
		}
		//if(cur_s<15)return;
		if(cur_s<10)return;
		dock = 1;
		sbc1 = rawline.substr(10, 7);
		umi = rawline.substr(0, 10);
		if(t==-1){
			umi = "N" + rawline.substr(0, 9);
		}
		else if(t==-2){
			umi = "NN" + rawline.substr(0, 8);
		}

		//2nd bc
		bait = "ATCCACGTGCTTGAG";
		cur_s = 0;
		int tt = t;
		for(int i = -2; i < 3; ++i){
			string qu = rawline.substr(54+i+tt, 15);
			int score = align_score(qu, bait);
			if(score < cur_s)continue;
			cur_s = score;
			tt = t + i;
		}
		//if(cur_s<13)return;
		if(cur_s<7)return;
		t = tt;
		sbc2 = rawline.substr(47+t, 7);
		dock = 2;

		//4th BC
		bait = "GCATTCGAG";
		cur_s = 0;
		for(int i = -2; i < 3; ++i ){
			string qu = rawline.substr(69+i+tt, 9);
			int score = align_score(qu, bait);
			if(score < cur_s)continue;
			cur_s = score;
			tt = t + i;
		}
		//if(cur_s<9)return;
		if(cur_s<5)return;
		t = tt;
		if(rawline.length() < 83+t)return;
		//sbc4 = rawline.substr(86+t, 3);
		sbc4 = rawline.substr(78+t, 5);
		//bsbc4 = rawline.substr(82+t, 7);
		dock = 4;		

		type= "a";
		return;
	}

	
	bool is_valid(){
		return valid;
	}
	int where_dock(){
		return dock;
	}
};


#endif //REACHTOOLS_CXSTRING_H
