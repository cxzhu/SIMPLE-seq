//
// Chenxu Zhu (chz272@ucsd.edu)
// 3/11/2020
//

#include <iostream>
#include "cxstring.h"
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <unistd.h>

using namespace std;


void help();
void combine(string sample);
void combine2(string sample, string ty);
void combine2_help();


string version = "2020.03.11";

int main(int argc, char** argv) {
	if(argc < 2){
		help();
		return 1;
	}
	string mod(argv[1]);
	if(mod == "-h" || mod == "help" || mod == "--help"){
		help();
		return 0;
	}

	if(mod == "combine"){
		if(argc < 3){
			combine2_help();
			return 1;
		}
		if(argc < 4){
			combine2(argv[2], "gz");
			return 0;
		}
		combine2(argv[2], argv[3]);
		return 0;
	}

	
	if(mod == "convert"){
		if(argc < 3){
			convert2_help();
			return 1;
		}
		convert2(argv[2]);
		return 0;
	}

	return 0;
}

//  local functions
void help(){
	cout << "Processiong of SIMPLE-seq data, modified from Paired-seq/Tag 2-round ligation barcodes." << "\nVersion: " << version << " (cxzhu@pku.edu.cn)" << endl;
	return;
}


// processing SIMPLE-seq with 2-round ligation
void combine2(string r2, string ty){
	int total = 0;
	int pass = 0;
	string s1;
	string s2;
	string s3;
	if(ty == "gz"){
		s1 = "zcat ";
		s2 = r2 + "_R1.fq.gz";
	}
	else if(ty == "bz2"){
		s1 = "bzcat ";
		s2 = r2 + "_R1.fastq.bz2";
	}
	s3 = s1 + s2;
	FILE * red1;
	red1 = popen(s3.c_str(), "r");
	s1 = "gzip - > ";
	s2 = r2 + "_combined.fq.gz";
	s3 = s1 + s2;
	FILE * outfile;
	outfile = popen(s3.c_str(), "w");
	if(ty == "gz"){
		s1 = "zcat ";
		s2 = r2 + "_R2.fq.gz";
	}
	else if(ty == "bz2"){
		s1 = "bzcat ";
		s2 = r2 + "_R2.fastq.bz2";
	}
	s3 = s1 + s2;
	FILE * red2;
	red2 = popen(s3.c_str(), "r");
	char buffer[2000];
	fqline in_line1;
	fqline in_line2;
	while(fgets(buffer, sizeof(buffer), red1)){
		++total;
		string line1(buffer);
		fgets(buffer, sizeof(buffer), red2);
		string line2(buffer);
		line1 = cxstring::chomp(line1);
		line2 = cxstring::chomp(line2);
		in_line1.read_part_record(red1, line1);
		in_line2.read_part_record(red2, line2);
		read2_2r read_2;
		read_2.init(in_line2.seq);
		read_2.trim();
		string new_seq = read_2.sbc1 + read_2.sbc2 + read_2.sbc4;
		string umi = read_2.umi;
		in_line2.seq = new_seq;
		in_line2.qual = in_line2.qual.substr(0, in_line2.seq.length());
		if(in_line2.seq.length()!=19)continue;
		string a;
		stringstream as;
		as << in_line2.readname;
		as >> a;
		in_line2.readname = a + ":" + umi + ":" + in_line1.seq + ":" + in_line1.qual;
		in_line2.write_record(outfile);
		++pass;
	}
	pclose(red1);
	pclose(red2);
	pclose(outfile);

	cout << endl << "Sample: " << r2 << endl << "Total reads: " << total << endl << "PF reads: " << pass << endl;
	return;

}
void combine2_help(){
	cout << "simpleconv combine sample_Prefix gz/bz2(default gz)" << endl;
}


