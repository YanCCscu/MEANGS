#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <algorithm>
#include <numeric>
#include <vector>
#include <unistd.h>
#include <time.h>
#include "pcre++.h"
using namespace std;
using namespace pcrepp;
using std::ifstream;
using std::ofstream;


string opt_f;
int opt_m;
int opt_o;
int opt_v;
int opt_r;
int opt_p;
int opt_d;
int opt_l;
int opt_a;
int opt_z;
int opt_e;
string opt_g;
string opt_s;
int opt_t;
string opt_b;
int opt_c;
int opt_x;
int opt_n;
int opt_h;
int opt_i;
int opt_w;
int opt_j;
int opt_q;
int opt_y;
int opt_u;

int targetwordlen = 15;
int base_overlap = 2;
int min_overlap = 20;
int verbose = 0;
int MIN_READ_LENGTH = 21;
int SEQ_SLIDE = 1;
double min_base_ratio = 0.7;
int paired = 0;
int min_links = 5;
double max_link_ratio = 0.3;
int contig_size_cutoff = 100;
double insert_stdev = 0.75;
string unpaired_file = "no-g";
string seed_file = "no-targetSeed";
int max_trim = 0;
string base_name = "";
int tracked = 0;
int forcetrack = 0;
int max_count_trim = 10;
int min_tig_overlap = 20;
int npad_gaps = 0;
int ignorehead = 0;
int space_restriction = 0;
int min_depth_of_coverage = 0;
int tie_breaker = 0;
int ignore_read = 0;
int independent = 1;


map<int, int> per;
int MAX = 0;
int MAX_TOP = 1500; //maxmum anchoring edges
int TRACK_COUNT = 0;
int LR_MAXDIST_APART = 60000;
int illuminaLengthCutoff = 500;
int LOW_COVERAGE_TIG_IN_A_ROW = 100; // -w set determinate value
string assemblyruninfo = "";

ofstream LOG;
ofstream SHO;
ofstream SC;


struct SEEDDATA
{
	int count;
	map<string, string> names;
	string seed_name;
	string ori;
};

unordered_map<string, SEEDDATA> seed;


struct MATEPAIRDATA
{
	int bt;
	int is;
};

struct SETDATA
{
	int count;
	map<string, string> names;
	int grace;
};

struct TRACKDATA
{
	int start;
	int end;
	map<string, string> names;
	int cov;
};

struct TRACKALLDATA
{
	int start;
	int end;
	int tig;
};

struct BARCODEDATA
{
	int lastlength;
	int count;
};

struct PSRDATA
{
	int start;
	int end;
};

struct PILEUPDATA
{
	string tig;
	string interest;
	string qua;
	string seq;
};

struct PAIRDATA
{
	int links;
	int gaps;
};


Pcre reg_globle;
bool regexMatch(const string &str,const string &expression,const string &flag)
{
	reg_globle = Pcre(expression,flag);
	if (reg_globle.search(str))
	{
		return true;
	}
	else
		return false;
}
bool regexMatch(const string &str, const string &expression)
{
	reg_globle = Pcre(expression);
	if (reg_globle.search(str))
	{
		return true;
	}
	else
		return false;
}

string  getCurrentTimeStr()
{
	struct tm t;  
	time_t now;  
	time(&now);  
	localtime_r(&now, &t); 
	stringstream ss;
	ss << t.tm_year + 1900 << "-" << t.tm_mon + 1 << "-" << t.tm_mday << "-" << t.tm_hour << ":" << t.tm_min << ":" << t.tm_sec;
	string rt = ss.str();
	ss.clear();
	return rt;
}

unordered_map<string, string> encodeBases()
{
	unordered_map<string, string> encoded;
	vector<string> pos1 = { "A","C","G","T" };
	vector<string> pos2 = pos1;
	vector<string> pos3 = pos1;
	vector<string> pos4 = pos1;
	string name1 = "的一国在人了有中是年和大业不为发会工经上地市要个产这出行作生家以成到日民来我部对进多全建他公开们场展时理新方主企资实学报制政济用同于法高长现本月定化加动合品重关机分力自外者区能设后就等体下万元社过前面农也得与说之员而务利电文事可种总改三各好金第司其从平代当天水省提商十管内小技位目起海所立已通入量子问度北保心还科委都术使明着次将增基名向门应里美由规今题记点计去强两些表系办教正条最达特革收二期并程厂如道际及西口京华任调性导组东路活广意比投决交统党南安此领结营项情解议义山先车然价放世间因共院步物界集把持无但城相书";
	vector<string> chararr;
	for (int i = 0; i < name1.size(); i += 2)
	{
		string tn = name1.substr(i, 2);
		int j;
		for (j = 0; j < chararr.size(); j++)
		{
			if (tn == chararr[j])
				break;
		}
		if (j == chararr.size())
		{
			chararr.push_back(tn);
		}
	}
	
	
	int el = 0;
	for (int p1 = 0; p1 < pos1.size(); p1++)
	{
		for (int p2 = 0; p2 < pos2.size(); p2++)
		{
			for (int p3 = 0; p3 < pos3.size(); p3++)
			{
				for (int p4 = 0; p4 < pos4.size(); p4++)
				{
					string quad = pos1[p1] + pos2[p2] + pos3[p3] + pos4[p4];
					encoded[quad] = chararr[el];
					el++;
				}
			}
		}
	}
	return encoded;

}

string reverseComplement(string word)
{
	transform(word.begin(), word.end(), word.begin(), ::toupper);
	reverse(word.begin(), word.end());
	for (int i = 0; i < word.length(); i++)
	{
		if (word[i] == 'A')
			word[i] = 'T';
		else if (word[i] == 'T')
			word[i] = 'A';
		else if (word[i] == 'G')
			word[i] = 'C';
		else if (word[i] == 'C')
			word[i] = 'G';
	}
	return word;
}

void loadSeed(string file, int targetwordlen, unordered_map<string, SEEDDATA> &seed, unordered_map<string, int> &seedsplit)
{
	ifstream IN;
	IN.open(file); 
	if (!IN)
	{
		perror(("Can't open " + file).c_str());
	}
	string subseq = "";
	string prev = "";
	while (!IN.eof())
	{
		string lineStr;
		getline(IN, lineStr);

		//Pcre reg("^\\>(\\S+)");
		//Pcre reg1("(\\S+)", "i");
		if (regexMatch(lineStr,"^\\>(\\S+)"))
		{
			string head;
			if (reg_globle.matches() >0)
			{
				head = reg_globle.get_match(0);
			}
			int subseq_length = subseq.length();
			if (subseq_length > MAX)
			{
				MAX = subseq_length;
			}
			if ((head != prev) && (subseq != "") && (subseq_length >= MIN_READ_LENGTH) && (subseq_length >= min_overlap))
			{
				string ucsub = "";
				ucsub.resize(subseq.size());
				transform(subseq.begin(), subseq.end(), ucsub.begin(), ::toupper);
				seed[ucsub].count++;
				seed[ucsub].names[prev] = "";
				seed[ucsub].seed_name = prev;
				seed[ucsub].ori = subseq;
				for (int pos = 0; pos <= subseq_length - targetwordlen; pos++)
				{
					string word = ucsub.substr(pos, targetwordlen);
					string word_rc = reverseComplement(word);
					seedsplit[word] = 1;
					seedsplit[word_rc] = 1;
				}

				Pcre reg_subseq("([BDEFHIJKLMNOPQRSUVWXYZ])","i");
				if (reg_subseq.search(subseq))
				{
					cout << "WARNING: sequence >" <<prev<<" contains non-AGCG characters"<< reg_subseq.get_match(0)<<".\n" << endl;
				}
			}
			subseq = "";
			prev = head;
		}
		else if (regexMatch(lineStr, "(\\S+)", "i"))
		{
			if (reg_globle.matches() >0)
			{
				subseq += reg_globle.get_match(0);
			}
			
		}
	}
	int subseq_length = subseq.length();
	if (subseq_length > MAX)
		MAX = subseq_length;
	if (subseq != "" && subseq_length >= MIN_READ_LENGTH && subseq_length >= min_overlap)
	{
		string ucsub = "";
		ucsub.resize(subseq.size());
		transform(subseq.begin(), subseq.end(), ucsub.begin(), ::toupper);
		seed[ucsub].count++;
		seed[ucsub].names[prev] = "";
		seed[ucsub].seed_name = prev;
		seed[ucsub].ori = subseq;
		for (int pos = 0; pos <= subseq_length - targetwordlen; pos++)
		{
			string word = ucsub.substr(pos, targetwordlen);
			string word_rc = reverseComplement(word);
			seedsplit[word] = 1;
			seedsplit[word_rc] = 1;
		}
		Pcre reg_subseq("([BDEFHIJKLMNOPQRSUVWXYZ])", "i");
		if (reg_subseq.search(subseq))
		{
			cout << "WARNING: the fasta sequence > "<<prev<<" in your seed file contains characters other than ACGTand may prevent proper contig extension."<< endl;
		}
	}
	IN.close();

}


void getSub(string str, vector<string> &vs)
{
	int t = str.length()/4;
	for (int i = 0; i < t; i++)
	{
		int idx = i * 4;
		vs.push_back(str.substr(idx, 4));
	}
}

void loadSequence(unordered_map<string, SETDATA> &set, unordered_map<string,map<string,string>> &bin, int &ctrd, string seq,  unordered_map<string, string> &e, string head, int up_iz, const unordered_map<string, int> &seedsplit, int space_restriction, int targetwordlen, int &recruit_mate)
{

	string orig="";
	orig.resize(seq.size());
	transform(seq.begin(), seq.end(), orig.begin(), ::toupper);
	int orig_mer = orig.length();

	if (orig != "" && orig_mer >= MIN_READ_LENGTH && orig_mer >= min_overlap)
	{
		if ((space_restriction) && (targetwordlen > orig_mer))
		{
			string msg = "ERROR loading reads:" + orig + "is too short" + to_string(orig_mer) + " as -j " + to_string(orig_mer) + "specified.\n";
			LOG << msg << endl;
			perror(msg.c_str());
		}

		vector<string> f;
		getSub(orig, f);
		string rc = reverseComplement(orig);
		vector<string> r;
		getSub(rc, r);
		string first_f = orig.substr(0, targetwordlen);
		string first_r = rc.substr(0, targetwordlen);

		if (orig_mer > MAX)
			MAX = orig_mer;
		

		if (!space_restriction || (space_restriction && (seedsplit.find(first_f)!=seedsplit.end()&&(seedsplit.at(first_f)) || seedsplit.find(first_r) != seedsplit.end() && seedsplit.at(first_r))) || space_restriction && recruit_mate)
		{
			set[orig].count++;
			set[orig].names[head] = "";
			set[orig].grace = up_iz;
			string barcode = "";
			if (regexMatch(head, "\\_(\\w+)"))
			{
				if (reg_globle.matches() >= 1)
				{
					barcode = reg_globle.get_match(0);
				}
			}

			string ef = e[f[0]] + e[f[1]] + e[f[2]] + e[f[3]];
			bin[ef][orig] = barcode;
			string er = e[r[0]] + e[r[1]] + e[r[2]] + e[r[3]];
			bin[er][rc] = barcode;

			recruit_mate++;
			ctrd++;
			//cout << "\r" << ctrd;
			if (ctrd%100 == 0) cout << "\r" << ctrd;
			cout << "\r";
		}

		f.clear();
		f.shrink_to_fit();
		r.clear();
		r.shrink_to_fit();
	}
	else if (orig != "")
	{
		if (orig_mer < MIN_READ_LENGTH)
		{
			SHO << seq << " Input sequence shorter than minimum read length allowed (" << orig_mer << "<" << MIN_READ_LENGTH <<"nt)"<< endl;
		}
		else if (orig_mer < min_overlap)
		{
			//print SHO "$seq\tInput sequence shorter than minimum overlap specified($orig_mer < -m $min_overlap)\n";
			SHO << seq << " Input sequence shorter than minimum overlap specified " << orig_mer << "< -m" << min_overlap << endl;
		}
	}
	if (MAX > MAX_TOP)
		MAX = MAX_TOP;
}

void readFasta(unordered_map<string, map<string, MATEPAIRDATA>> &matepair, unordered_map<string, SETDATA> &set, unordered_map<string, map<string, string>> &bin, string file, string shortfile, int paired,  unordered_map<string, string> &encoded, const unordered_map<string, int> &seedsplit, int space_restriction, int targetwordlen, int ignorehead, int &sumall, int &ctall)
{
	int ctrd = 0;
	int ctline = 0;
	string head = "";
	int insert_size = 0;
	int up_iz = 0;

	ifstream IN;
	IN.open(file); 
	if (!IN)
	{
		perror(("Can't open " + file).c_str());
	}
	
	SHO.open(shortfile); 
	if (!SHO)
	{
		perror(("Can't open " + shortfile).c_str());
	}
	//cout << "Sequence reads loaded:" << endl;
	while (!IN.eof())
	{
		string lineStr;
		getline(IN, lineStr);
		ctline++;
		Pcre reg("^([^\\>]*)$","i");
		if (reg.search(lineStr))
		{
			string sdna = reg.get_match(0);
			int recruit_mate = 0; // to track a read matches a target
			if (paired)
			{
				Pcre reg_sdna("([ACGT]*)\\:([ACGT]*)","i");
				if (reg_sdna.search(sdna))
				{
					// don't crash on Ns, but ignore see below
					string rd1 = reg_sdna.get_match(0);
					string rd2 = reg_sdna.get_match(1);
					string head1 = head + "1";
					string head2 = head + "2";
					//#num
					int len = rd1.length() + rd2.length();
					sumall += len;
					ctall += 2;

					if (ignorehead)
						head1 = "";
					if (ignorehead)
						head2 = "";
					loadSequence(set, bin, ctrd, rd1, encoded, head1, up_iz, seedsplit, space_restriction, targetwordlen, recruit_mate);
					loadSequence(set, bin, ctrd, rd2, encoded, head2, up_iz, seedsplit, space_restriction, targetwordlen, recruit_mate);
					if (recruit_mate == 1)
						loadSequence(set, bin, ctrd, rd1, encoded, head1, up_iz, seedsplit, space_restriction, targetwordlen, recruit_mate); // re - visit the first read if second has match in target space
					if (recruit_mate == 2)
					{
						matepair[rd1][rd2].is = insert_size;
						matepair[rd2][rd1].is = insert_size;
						matepair[rd1][rd2].bt = 0;
						matepair[rd2][rd1].bt = 0;
					}
				}
				else
				{
					string pairing_failure_message = "Input error at line "+to_string(ctline)+": The sequence \""+sdna+"\" is not in the right format for paired-end reads  -- Fatal\nMake sure your input is in the form (input sequences can be of variable lengths):\n\n>test\nGCTACGACTATGACATACAGT:GTAGATTGATCGCATGCACGCT\n\nWhere : separates paired reads.  Spaces, <<.>> or any characters other than A,C,G or T in your input file might have caused this error, including reads with Ns.\n";
					cout << pairing_failure_message;
					LOG << pairing_failure_message << endl;
					LOG.close();
					assemblyruninfo += pairing_failure_message + "\n";
					exit(0);
				}
			}
			else
			{
				if (ignorehead)
					head = "";
				if (regexMatch(sdna, "^([ACGT]*)$", "i"))
					loadSequence(set, bin, ctrd, sdna, encoded, head, up_iz, seedsplit, space_restriction, targetwordlen, recruit_mate);
				int len = sdna.length();
				sumall += len;
				ctall++;
			}

		}
		else if (regexMatch(lineStr, "^\\>(\\S+)"))
		{
			if (reg_globle.matches() >= 1)
			{
				head = reg_globle.get_match(0);
			}
			if (paired) 
			{
				if (regexMatch(head, "(\\S+)\\:(\\d+)$"))
				{
					if (reg_globle.matches() >= 2)
					{
						head = reg_globle.get_match(0);
						insert_size = stoi(reg_globle.get_match(1));
					}
				}
				int min_allowed = -1 * (insert_stdev * insert_size);
				up_iz = insert_size - min_allowed;
			}

			if (head == "" || (paired && insert_size == 0))
			{
				string input_failure_message = "Input error at line #" + to_string(ctline) + ":"+ lineStr+" -- Either you forgot to name the mate pair template (head=" + head + "), improperly formatted the insert size (iz=" + to_string(insert_size) + ") # or both.  Headers should look like: >templateA:200\nIf you are not using the -p 1 option, you can leave \":insert_size\" out.\n";
				cout << input_failure_message;
				LOG << input_failure_message << endl;
				assemblyruninfo += input_failure_message + "\n";
				LOG.close();
				exit(0);
			}
		}

	}

	IN.close();
	SHO.close();
	string read_number_message = "total sequences: "+to_string(ctrd)+"\tunique: "+to_string(set.size())+"\n";
	//cout<< read_number_message;
	LOG << read_number_message << endl;
	assemblyruninfo+= read_number_message + "\n";


}




void deleteData(unordered_map<string, map<string, string>>  &bin, unordered_map<string, SETDATA> &set, unordered_map<string, SEEDDATA> &seed, string sequence,  unordered_map<string, string> &encoded)
{
	vector<string> o;
	getSub(sequence, o);
	string comp_seq = reverseComplement(sequence);
	vector<string> c;
	getSub(comp_seq, c);
	//remove k - mer from hash table and prefix tree
	string eo = encoded[o[0]] + encoded[o[1]] + encoded[o[2]] + encoded[o[3]];
	string ec = encoded[c[0]] + encoded[c[1]] + encoded[c[2]] + encoded[c[3]];
	bin[eo].erase(sequence);
	bin[ec].erase(comp_seq);

	set.erase(sequence);
	seed.erase(sequence);
	

}



vector<string> split(const string &s, const string &seperator)
{
	vector<string> result;
	typedef string::size_type string_size;
	string_size i = 0;

	while (i != s.size())
	{
		//charactor not seperator
		int flag = 0;
		while (i != s.size() && flag == 0)
		{
			flag = 1;
			for (string_size x = 0; x < seperator.size(); ++x)
				if (s[i] == seperator[x])
				{
					++i;
					flag = 0;
					break;
				}
		}

		//extract contents between seperators
		flag = 0;
		string_size j = i;
		while (j != s.size() && flag == 0)
		{
			for (string_size x = 0; x < seperator.size(); ++x)
				if (s[j] == seperator[x])
				{
					flag = 1;
					break;
				}
			if (flag == 0)
				++j;
		}
		if (i != j)
		{
			result.push_back(s.substr(i, j - i));
			i = j;
		}
	}
	return result;
}

void collectOverhang(unordered_map<int, map<string, int>> &overhang, string newpass, string dangle, unordered_map<string, SETDATA> &set, int verbose)
{

	string over = dangle;
	int ct_oh = 0;

	//foreach string bz(@over)
	for (int i = 0; i < over.size(); i++)
	{
		string bz;
		bz.push_back(over[i]);
		ct_oh++; //tracks overhang position passed the seed
		overhang[ct_oh][bz] += set[newpass].count;   //reflects read coverage(often real duplicates)
		if (verbose)
			cout<< ct_oh<< "-"<< bz <<"="<< overhang[ct_oh][bz]<<"\n";
	}
}



void doExtension(string direction, int orig_mer, string &seq, unordered_map<string, SETDATA> &set, unordered_map<string, map<string, string>> &bin, int &reads_needed, int &total_bases, int min_overlap, int base_overlap, double min_base_ratio, int verbose, unordered_map<string, TRACKDATA> &track, int paired, int tig_count, int max_trim, unordered_map<string, string> encoded, const unordered_map<string, map<string, MATEPAIRDATA>> &matepair, int tie_breaker, int ignore_read, unordered_map<string, BARCODEDATA> &barcodeList)
{
	int extended = 1;
	int trim_ct = 0;     //trim counter - keeps track of 3'-end trim
	if (orig_mer > MAX)
	{
		orig_mer = MAX;
	} //Deals with special cases where the seed sequences are different from the read set(and possibly very large) - goal here is not to increase sequence coverage of seed, but rather to extend it.

	//#TRIM:
	while (trim_ct <= max_trim)
	{
		while (extended)
		{
			int growing_tig_length = seq.length();
			int pos = 0;
			int span = 0;
			if (growing_tig_length >= MAX)
			{
				//seq is length of contig being extended-- if larger than largest read, make sure the largest read could align and all subsequent rds.
				span = MAX - TRACK_COUNT;
			}
			else
			{
				span = growing_tig_length - TRACK_COUNT;
			}

			while (span >= min_overlap)
			{
				// will slide the subseq, until the user - defined min overlap size
				pos = growing_tig_length - span;
				if (verbose)
					LOG << "MAX:" << MAX << ",SPAN:" << span << ", POS:" << pos << endl;

				string subseq = seq.substr(pos, span);              //make a sub - sequence of length l - (1..i) for searching
				vector<string> s;
				getSub(subseq, s);
				string es = encoded[s[0]] + encoded[s[1]] + encoded[s[2]] + encoded[s[3]];
				
				map<string, string> subset = bin[es]; //Will grab everything even the reverse complement ones
				if (verbose)
					cout << direction << "SEARCH Position:" << pos << " Span:" << span << " - Subseq:" << subseq << "Previous:" << seq << endl;

				vector<std::pair<string, string>> tmpSubset;
				for (auto& siter : subset)
					tmpSubset.push_back(siter);
				sort(tmpSubset.begin(), tmpSubset.end(),
					[=](std::pair<string, string>& a, std::pair<string, string>& b) { return a.second > b.second; });

				// SEARCH -- this cycles through limited k-mer space
				for (int i = 0; i < tmpSubset.size(); i++)
				{
					string pass = tmpSubset[i].first;
					string ex1 = "^" + subseq + "([ACGT]+)";
					string ex2 = pass;
					if(regexMatch(pass,ex1)||regexMatch(subseq,ex2))
					{
						// OVERHANG || EMBEDDED
						//print "BC subset->{pass}\n";
						string sp = subset[pass];
						if (barcodeList.find(sp) != barcodeList.end())
						{
							// barcode seen before
							if ((seq.length() - barcodeList[sp].lastlength) > LR_MAXDIST_APART)
							{
								barcodeList[sp].count = 0; ////// reset barcode to zero, too far apart, new molecule
							}
						}
						barcodeList[sp].lastlength = seq.length(); ////// new tracks all barcodes encountered
						barcodeList[sp].count++; ////// new tracks all barcodes encountered
					}
				}
				span--;
				
				//map<string, string> dsubset;
				//subset.swap(dsubset);

			}//while overlap >= user - defined - m minimum
			// END NEW ADDITION==========================================
			pos = 0;
			span = 0;
			if (growing_tig_length >= MAX)
			{
				// seq is length of contig being extended-- if larger than largest read, make sure the largest read could align and all subsequent rds.
				span = MAX - TRACK_COUNT;
			}
			else
			{
				span = growing_tig_length - TRACK_COUNT;
			}

			unordered_map<int, map<string, int> >   overhang;
			map<string, int> overlapping_reads;
			int mlong = 0;

			for (int x = 1; x <= (orig_mer * 2); x++)
			{
				map<string, int> omap;
				omap["A"] = 0;
				omap["C"] = 0;
				omap["G"] = 0;
				omap["T"] = 0;
				overhang[x] = omap;
			}


			//COLLECT SEQUENCES 
			while (span >= min_overlap)
			{
				// will slide the subseq, until the user - defined min overlap size
				pos = growing_tig_length - span;
				if (verbose)
					cout << "MAX:" << MAX << ",SPAN:" << span << ", POS:" << pos << endl;
				string subseq = seq.substr(pos, span);              //make a sub - sequence of length l - (1..i) for searching
				vector<string> s;
				getSub(subseq, s);
				string es = encoded[s[0]] + encoded[s[1]] + encoded[s[2]] + encoded[s[3]];

				map<string, string> subset = bin[es]; //Will grab everything even the reverse complement ones
				if (verbose)
					cout << direction << "SEARCH Position:" << pos << " Span:" << span << " - Subseq:" << subseq << "Previous:" << seq << endl;

				vector<std::pair<string, string>> tmpSubset;
				for (auto& siter : subset)
					tmpSubset.push_back(siter);
				sort(tmpSubset.begin(), tmpSubset.end(),
					[=](std::pair<string, string>& a, std::pair<string, string>& b) { return a.second > b.second; });
				
				// SEARCH -- this cycles through limited k-mer space
				for (int i = 0; i < tmpSubset.size(); i++)
				{
					string pass = tmpSubset[i].first;
					
					string ex1 = "^" + subseq + "([ACGT]+)";
					if (regexMatch(pass, ex1))
					{
						// OVERHANG 
						//can we align perfectly that subseq to another rd start ?
						string dangle;
						if (reg_globle.matches() > 0)
						{
							dangle = reg_globle.get_match(0);
						}
						//if(verbose)
							//cout<<"="<<
						unordered_map<string, PSRDATA> psr;
						string pass_rc = reverseComplement(pass);

						if (set.find(pass) != set.end())
						{
							psr[pass].start = pos + 1;
							psr[pass].end = pos + pass.length();
						}

						if (set.find(pass_rc) != set.end())
						{
							psr[pass_rc].start = pos + pass_rc.length();
							psr[pass_rc].end = pos + 1;
						}

						
						// CONSIDER CERTAIN READS FOR OVERLAP, PREFERABLY THOSE WITH LOGICAL MATES AND FWD READS IN TIG LARGE ENOUGH TO HAVE SUCH PAIRS
						for (auto &iter : psr)
						{
							string newpass = iter.first;
							if ((seq.length() <= 2000 || ((seq.length() > 2000 || seq.length() <= 10000) && barcodeList[subset[newpass]].count > 4) || (seq.length() > 10000 && barcodeList[subset[newpass]].count >= 9)) && paired && matepair.find(newpass)!=matepair.end() && (iter.second.end < iter.second.start) && (iter.second.start >= set[newpass].grace))
							{
								//paired, pairingRds, < -- - , outside grace
								map<string, MATEPAIRDATA> mateshash;
								if (matepair.find(newpass) != matepair.end())
								{
									mateshash = matepair.at(newpass);
								}
								
								//#MATESEARCH:
								for (auto &mIter : mateshash)
								{
									string matchingread = mIter.first;
									if (track.find(matchingread) != track.end())
									{
										//a mate has been found on this contig
										int insert_size = mIter.second.is;
										int min_allowed = -1 * (insert_stdev * insert_size);
										int low_iz = insert_size + min_allowed;
										int up_iz = insert_size - min_allowed;
										int A_start = track[matchingread].start;
										int A_end = track[matchingread].end;
										int pet_size = iter.second.start - A_start;
										//print "\tnewpass (psr->{newpass}{'start'} - psr->{newpass}{'end'}) and matchingread (A_start - A_end ) are on same tig#tig_count  [low_iz to up_iz]\n" if (verbose);
										if ((iter.second.start > A_start) && (A_start < A_end) && (pet_size >= low_iz) && (pet_size <= up_iz))
										{
											collectOverhang(overhang, newpass, dangle, set, verbose);
											overlapping_reads[pass]++;
											if (newpass.length() > illuminaLengthCutoff)
												mlong = 1;
											break;
										}
									}
								}
								//map<string, MATEPAIRDATA> dmateshash;
								//mateshash.swap(dmateshash);
							}
							else if (seq.length() <= 2000 || ((seq.length() > 2000 || seq.length() <= 10000) && barcodeList[subset[newpass]].count > 4) || (seq.length() > 10000 && barcodeList[subset[newpass]].count >= 9))
							{
								// not paired or paired (B-->): collect all overlaps
								collectOverhang(overhang, newpass, dangle, set, verbose);
								overlapping_reads[pass]++;
								if (newpass.length() > illuminaLengthCutoff)
									mlong = 1;
							}
						}//for $newpass

						//unordered_map<string, PSRDATA> dpsr;
						//psr.swap(dpsr);

					}
					else if (regexMatch(subseq,pass))
					{
						// EMBEDDED
						string complement_pass = reverseComplement(pass);
						
						//print "$pass found in $subseq ($set->{$pass}{'count'}) - deleting read: $pass and complement ($set->{$complement_pass}): $complement_pass\n\n" if ($verbose);
						if (set.find(pass) != set.end())
						{
							int current_reads = set[pass].count;
							int current_bases = pass.length() * current_reads;
							reads_needed += current_reads;
							total_bases += current_bases;
							if (paired)
							{
								track[pass].start = pos + 1;
								track[pass].end = pos + pass.length();
								track[pass].cov = current_reads;
								track[pass].names = set[pass].names;
							}
							deleteData(bin, set, seed, pass, encoded);
							//print "EMBED .. pass (current_reads)\n";
						}
						if (set.find(complement_pass) != set.end())
						{
							int current_reads = set[complement_pass].count;
							int current_bases = complement_pass.length() * current_reads;
							reads_needed += current_reads;
							total_bases += current_bases;
							if (paired)
							{
								track[complement_pass].end = pos + 1;
								track[complement_pass].start = pos + complement_pass.length();
								track[complement_pass].cov = current_reads;
								track[complement_pass].names = set[complement_pass].names;
							}
							deleteData(bin, set, seed, complement_pass, encoded);
							//print "EMBED .. complement_pass (current_reads)\n";
						}
					}
				}
				span--;
				//map<string, string> dsubset;
				//subset.swap(dsubset);
			}//while overlap >= user-defined -m minimum

			
			string consensus = "";
			if (verbose)
				cout << "Finished Collecting Overlapping Reads - BUILDING CONSENSUS...\n";
			//print Dumper(%$overlapping_reads) if ($verbose);
			int tmp_base_overlap = base_overlap;
			if (mlong)
			{
				tmp_base_overlap = 2;
			}

			//Build consensus
			vector<int> keys;
			for (auto const &iter : overhang)
			{
				int key = iter.first;
				keys.push_back(key);
			}
			sort(keys.begin(), keys.end());
			//#CONSENSUS
			for (int oi = 0; oi < keys.size(); oi++)
			{
				bool flag = true;
				int ohpos = keys[oi];
				//LOG << ohpos <<"A:"<< overhang[ohpos]["A"] << " C:" << overhang[ohpos]["C"] << " G:" << overhang[ohpos]["G"] << " T:" << overhang[ohpos]["T"] <<endl;
				if (ohpos)
				{

					int coverage = overhang[ohpos]["A"] + overhang[ohpos]["C"] + overhang[ohpos]["G"] + overhang[ohpos]["T"];
					//print "pos:$ohpos cov:$coverage A:$overhang->{$ohpos}{'A'} C:$overhang->{$ohpos}{'C'} G:$overhang->{$ohpos}{'G'} T:$overhang->{$ohpos}{'T'}\n" if ($verbose);
					if (coverage < tmp_base_overlap)
					{
						//print "COVERAGE BELOW THRESHOLD: $coverage < -o $tmp_base_overlap @ $ohpos :: will extend by: $consensus\n" if ($verbose);
						break;
					}
					map<string, int> baselist = overhang[ohpos];
					int ct_dna = 0;
					string previous_bz = "";


					vector<std::pair<string, int>> tempbl;
					for (auto& biter : baselist)
						tempbl.push_back(biter);
					sort(tempbl.begin(), tempbl.end(),
						[=](std::pair<string, int>& a, std::pair<string, int>& b) { return a.second > b.second; });

					//#BASE
					for (int bi = 0; bi < tempbl.size(); bi++)
					{
						string bz = tempbl[bi].first;
						//print "\t$ct_dna -> $bz..$baselist->{$previous_bz} > $baselist->{$bz}\n";
						if (ct_dna)
						{
							// the two most abundant bases at that position
							//print "\t\t$ct_dna\n";
							if (previous_bz != "" && ((double)baselist[previous_bz]/(double)coverage) >= min_base_ratio && baselist[previous_bz] > baselist[bz])
							{
								//# a simple consensus btw top 2 
								consensus += previous_bz;    //# build consensus
								//print "Added base $previous_bz (cov = $baselist->{$previous_bz}) to $consensus **\n" if ($verbose);
								break;
							}
							else
							{
								if (previous_bz != "" && tie_breaker && coverage <= 2)
								{
									//# testing limit coverage
									consensus += previous_bz;                //# build consensus
									//print "Added base $previous_bz (cov = $baselist->{$previous_bz}) to $consensus (Forced by -q $tie_breaker) **\n" if ($verbose);
									break;
								}
								else
								{
									//print "ISSUES EXTENDING: best base = $previous_bz (cov=$baselist->{$previous_bz}) at $ohpos.  Second-Best: $bz (cov=$baselist->{$bz}) (ratio best=$baselist->{$previous_bz} / total=$coverage) >= $min_base_ratio (-r) -- will terminate with $consensus\n" if ($verbose);
									flag = false;
									break;
								}
							}
						}
						previous_bz = bz;
						ct_dna++;
					}
					
					//map<string, int> dbaselist;
					//baselist.swap(dbaselist);
					if (!flag)
						break;
				}
			}
			mlong = 0;

			if (ignore_read == 0)
			{
				//#new option to ignore read mapping
				//# deal with sequence reads making up the consensus/newly formed contig
				if (consensus != "")
				{
					//print "Will extend $seq\nwith: $consensus\n\n" if ($verbose);
					string temp_sequence = seq + consensus;
					int position_buffer = 0;
					if (growing_tig_length > MAX)
					{
						temp_sequence = seq.substr(growing_tig_length - MAX, MAX) + consensus;
						position_buffer = growing_tig_length - MAX;
					}
					int integral = 0;
					string temp_sequence_portion = "";
					for (int i = 0; i < temp_sequence.length(); i++)
					{
						temp_sequence_portion += "-";
					}

					//LOG << "overlapping:-----------" << endl;
					for (auto &oIter : overlapping_reads)
					{
						string ro = oIter.first;
						//LOG << "ro:" << ro << endl;
						int or_pos = -99;
						// want the last position
						if (regexMatch(temp_sequence, ro))
						{
							or_pos = temp_sequence.rfind(ro) + ro.size() - 1;
						}
						
						if (or_pos > 0)
						{
							//TRACK COVERAGE TO PREVENT FAULTY EXTENSIONS
							/*
							string linestring = "";
							for (int i = 0; i < ro.length(); i++)
							{
								linestring += ".";
							}
							temp_sequence_portion.replace(or_pos, ro.length(), linestring);
							*/
							int repos= or_pos + ro.size();
							int tempsize = temp_sequence_portion.size();
							for (int i = or_pos; i <= repos; i++)
							{
								if (i <= tempsize - 1)
								{
									temp_sequence_portion[i] = '.';
								}
								else
								{
									temp_sequence_portion.push_back('.');
								}
							}

							or_pos += 1;
							or_pos += position_buffer;
							string complement_ro = reverseComplement(ro);
							//LOG << "comro:" << complement_ro << endl;
							//print "$ro found in $seq ($set->{$ro}{'count'}) - deleting read: $ro and complement ($set->{$complement_ro}{'count'}): $complement_ro\n\n" if ($verbose);
							if (set.find(ro) != set.end())
							{
								int current_reads = set[ro].count;
								//#print "fwd SET:current_reads BIN subset[ro}\n";
								int current_bases = ro.length() * current_reads;
								integral += current_reads;
								reads_needed += current_reads;
								total_bases += current_bases;
								if (paired)
								{
									track[ro].start = or_pos - ro.length() + 1;
									track[ro].end = or_pos;
									track[ro].cov = current_reads;
									track[ro].names = set[ro].names;
									//#my lro = length(ro);
									//#print "OVER FWD\nseq (growing_tig_length)\ntemp_sequence\nconsensus\nro (current_reads)\tend=or_pos (track[ro}{'start'}-track[ro}{'end'})\n\n";
								}
								deleteData(bin, set, seed, ro, encoded);
							}

							if (set.find(complement_ro) != set.end())
							{
								int current_reads = set[complement_ro].count;
								//#print "rc SET:current_reads BIN subset_rc[complement_ro}\n";
								int current_bases = complement_ro.length() * current_reads;
								integral += current_reads;
								reads_needed += current_reads;
								total_bases += current_bases;
								if (paired)
								{
									track[complement_ro].end = or_pos - ro.length() + 1;
									track[complement_ro].start = or_pos;
									track[complement_ro].cov = current_reads;
									track[complement_ro].names = set[complement_ro].names;
									//my $lro = length($ro);
									//print "OVER REV\n$seq ($growing_tig_length)\n$temp_sequence\n$consensus\n$complement_ro ($current_reads)\tstart=$or_pos($track->{$complement_ro}{'start'}-$track->{$complement_ro}{'end'})\n\n";
								}
								deleteData(bin, set, seed, complement_ro, encoded);
							}
						}
					}
					if (integral < tmp_base_overlap)
					{
						//print "No overlapping reads agree with the consensus or number of agreeing reads is lower than target coverage (tmp:tmp_base_overlap  -o base_overlap). Stopping extension" if (verbose);
						extended = 0;
					}
					else
					{
						if (regexMatch(temp_sequence_portion,"\\.(\\-*)?$"))
						{
							// will not extend a contig with 3' consensus bases if no reads overlap them (mitigate assembly errors)
							string str;
							if (reg_globle.matches() > 0)
							{
								str = reg_globle.get_match(0);
							}
							consensus = consensus.substr(0, consensus.length() - str.length());
						}
						seq += consensus;   // contig extension
						if (verbose)
							cout << "New Contig is: " << seq << "\n";
						extended = 1;
					}
				}
				else
				{
					// no consensus built, will stop the extension
					extended = 0;
				}
			}
			else
			{
				//New option Dec2014 to ignore reads making up consensus
				if (consensus != "")
				{
					seq += consensus;   // contig extension
					if (verbose)
						cout << "ALLA New Contig is: " << seq << "\n";
					extended = 1;
				}
				else
				{
					extended = 0;
				}
			}
			if (verbose)
				cout << "IGNORE:" << ignore_read << "\n";

			//unordered_map<int, map<string, int> >   doverhang;
			//overhang.swap(doverhang);	
			//map<string, int> doverlapping_reads;
			//overlapping_reads.swap(doverlapping_reads);

		}//while get the OK for extension

		trim_ct++;
		if (trim_ct <= max_trim)
		{
			if (seq.length() <= MIN_READ_LENGTH)
				break; //#terminate assembly if trimming becomes too agressive
			seq = seq.substr(0, seq.size()-1);
			extended = 1;
			//print "\ndirection EXTENSION ROUND trim_ct COMPLETE UNTIL max_trim nt TRIMMED OFF => TRIMMED SEQUENCE:seq\n\n" if (verbose);
		}


	}// while trimming within bounds
	 // Adjust the position if tracking paired reads in assembly
	if (paired)
	{
		for (auto &tIter : track)
		{
			string rd = tIter.first;
			track[rd].start = seq.length() - tIter.second.start + 1;
			track[rd].end = seq.length() - tIter.second.end + 1;
		}
	}
	//print "\n*** NOTHING ELSE TO BE DONE IN $direction - PERHAPS YOU COULD DECREASE THE MINIMUM OVERLAP -m (currently set to -m $min_overlap) ***\n\n" if ($verbose);

}


void trackReads( unordered_map<string, TRACKDATA> &track, unordered_map<string, TRACKALLDATA> &track_all, unordered_map<int, map<int, map<string, int>>>& alternate, int tig_count)
{

	unordered_map<string, TRACKDATA>::iterator tIter;
	for (tIter=track.begin();tIter!=track.end();)
	{
		string rd = tIter->first;

		if (track_all.find(rd) == track_all.end())
		{
			track_all[rd].tig = tig_count;
			track_all[rd].start = tIter->second.start;
			track_all[rd].end = tIter->second.end;
			alternate[tig_count][tIter->second.start][rd]= tIter->second.end;
			track.erase(tIter++);
		}
		else
		{
			tIter++;
		}
	}


}

int getDistance(int insert_size, int length_i, int start_i, int start_j)
{
	int insert_span = (length_i - start_i) + start_j;
	int gap_or_overlap = insert_size - insert_span;

	return gap_or_overlap;
}

void pairContigs(unordered_map<string, map<string, PAIRDATA>> &pair,  unordered_map<string, map<string, MATEPAIRDATA>> &matepair, unordered_map<string, TRACKALLDATA> &track_all, unordered_map<int, int> &tig_length, string issues, string distribution, int verbose)
{
	int ct_illogical = 0;
	int ct_ok_contig = 0;
	int ct_ok_pairs = 0;
	int ct_problem_pairs = 0;
	int ct_iz_issues = 0;
	int ct_single = 0;
	int ct_both = 0;

	map<int, int>ct_both_hash;
	map<int, int>ct_ok_pairs_hash;
	map<int, int>ct_illogical_hash;
	map<int, int>ct_ok_contig_hash;
	map<int, int>ct_problem_pairs_hash;
	map<int, int>ct_iz_issues_hash;
	map<int, int>ct_single_hash;
	map<string, PAIRDATA>err;
	map<int, int>track_insert;

	if (verbose)
		cout << "Pairing contigs...\n";
	ofstream PET;
	PET.open(issues); 
	if (!PET)
	{
		PET.close();
		perror(("Can't write to " + issues).c_str());
	}
	for (auto &mIter : matepair)
	{
		string read_a = mIter.first;
		map<string, MATEPAIRDATA> mateslist = matepair[read_a];
		for (auto &lIter : mateslist)
		{
			string read_b = lIter.first;
			if (!matepair[read_a][read_b].bt && !matepair[read_b][read_a].bt)
			{
				matepair[read_a][read_b].bt = 1;
				matepair[read_b][read_a].bt = 1;
				int insert_size = mateslist[read_b].is;
				int min_allowed = -1 * (insert_stdev * insert_size);
				int low_iz = insert_size + min_allowed;
				int up_iz = insert_size - min_allowed;
				if (verbose)
					cout << "Pair read1=" << read_a << " read2=" << read_b << endl;
				if (track_all[read_a].tig && track_all[read_b].tig)  // both pairs assembled
				{
					ct_both++;
					ct_both_hash[insert_size]++;

					int tig_a = track_all[read_a].tig;
					int tig_b = track_all[read_b].tig;

					string ftig_a = "f" + to_string(tig_a);
					string ftig_b = "f" + to_string(tig_b);

					string rtig_a = "r" + to_string(tig_a);
					string rtig_b = "r" + to_string(tig_b);

					int A_length = tig_length[tig_a];
					int A_start = track_all[read_a].start;
					int A_end = track_all[read_a].end;

					int B_length = tig_length[tig_b];
					int B_start = track_all[read_b].start;
					int B_end = track_all[read_b].end;

					if (tig_a != tig_b)//paired reads located on <> contigs
					{
						//Determine most likely possibility
						if (A_start < A_end)
						{
							if (B_end < B_start)
							{
								int d = getDistance(insert_size, A_length, A_start, B_start);
								if (verbose)
									cout << "A-> <-B  WITH " << tig_a << " -> <- " << tig_b << "GAP " << d << "A=" << A_length << "(" << A_start - A_end << ")" << "B=" << B_length << "(" << B_start - B_end << ") Alen, Astart,Bstart\n";
								if (d >= min_allowed)
								{
									pair[ftig_a][ftig_b].links++;
									pair[ftig_a][ftig_b].gaps += d;
									pair[rtig_b][rtig_a].links++;
									pair[rtig_b][rtig_a].gaps += d;
									ct_ok_pairs++;
									ct_ok_pairs_hash[insert_size]++;
								}
								else
								{
									string err_pair = ftig_a + "-" + ftig_b;
									err[err_pair].links++;
									err[err_pair].gaps += d;
									ct_problem_pairs++;
									ct_problem_pairs_hash[insert_size]++;
									PET << "Pairs unsatisfied in distance within a contig pair.  A-> <-B  WITH tig#" << tig_a << "->" << d << "<- tig#r." << tig_b << ", A=" << A_length << "  nt (start:" << A_start << ", end:" << A_end << ") B=" << B_length << " nt (start:" << B_start << ", end:" << B_end << ") CALCULATED DISTANCE APART: " << d << "<" << min_allowed << endl;
								}


							}
							else
							{
								int rB_start = B_length - B_start;
								int d = getDistance(insert_size, A_length, A_start, rB_start);
								if (verbose)
									cout << "A-> <-rB  WITH " << tig_a << " -> <- r." << tig_b << "GAP " << d << "A=" << A_length << "(" << A_start - A_end << ")" << "B=" << B_length << "(" << B_start - B_end << ") Alen, Astart,rBstart\n";
								if (d >= min_allowed)
								{
									pair[ftig_a][rtig_b].links++;
									pair[ftig_a][rtig_b].gaps += d;
									pair[ftig_b][rtig_a].links++;
									pair[ftig_b][rtig_a].gaps += d;
									ct_ok_pairs++;
									ct_ok_pairs_hash[insert_size]++;
								}
								else
								{
									string err_pair = ftig_a + "-" + rtig_b;
									err[err_pair].links++;
									err[err_pair].gaps += d;
									ct_problem_pairs++;
									ct_problem_pairs_hash[insert_size]++;
									PET << "Pairs unsatisfied in distance within a contig pair.  A-> <-rB  WITH tig#" << tig_a << "->" << d << "<- tig#r." << tig_b << ", A=" << A_length << "  nt (start:" << A_start << ", end:" << A_end << ") B=" << B_length << " nt (start:" << B_start << ", end:" << B_end << ") CALCULATED DISTANCE APART: " << d << "<" << min_allowed << endl;
								}
							}
						}
						else
						{
							if (B_end > B_start)
							{
								int d = getDistance(insert_size, B_length, B_start, A_start);
								if (verbose)
									cout << "B-> <-A  WITH " << tig_b << " -> <- " << tig_a << "GAP " << d << "A=" << A_length << "(" << A_start - A_end << ")" << "B=" << B_length << "(" << B_start - B_end << ") Alen, Astart,Bstart\n";
								if (d >= min_allowed)
								{
									pair[ftig_b][ftig_a].links++;
									pair[ftig_b][ftig_a].gaps += d;
									pair[rtig_a][rtig_b].links++;
									pair[rtig_a][rtig_b].gaps += d;
									ct_ok_pairs++;
									ct_ok_pairs_hash[insert_size]++;
								}
								else
								{
									string err_pair = ftig_b + "-" + ftig_a;
									err[err_pair].links++;
									err[err_pair].gaps += d;
									ct_problem_pairs++;
									ct_problem_pairs_hash[insert_size]++;
									PET << "Pairs unsatisfied in distance within a contig pair.  B-> <-A  WITH tig#" << tig_b << "->" << d << "<- tig#r." << tig_a << ", B=" << B_length << "  nt (start:" << B_start << ", end:" << B_end << ") A=" << A_length << " nt (start:" << A_start << ", end:" << A_end << ") CALCULATED DISTANCE APART: " << d << "<" << min_allowed << endl;
								}


							}
							else
							{
								int rB_start = B_length - B_start;
								int d = getDistance(insert_size, B_length, rB_start, A_start);
								if (verbose)
									cout << "rB-> <-A  WITH " << tig_b << " -> <- r." << tig_a << "GAP " << d << "A=" << A_length << "(" << A_start - A_end << ")" << "B=" << B_length << "(" << B_start - B_end << ") Alen, Astart,rBstart\n";
								if (d >= min_allowed)
								{
									pair[rtig_b][ftig_a].links++;
									pair[rtig_b][ftig_a].gaps += d;
									pair[rtig_a][ftig_b].links++;
									pair[rtig_a][ftig_b].gaps += d;
									ct_ok_pairs++;
									ct_ok_pairs_hash[insert_size]++;
								}
								else
								{
									string err_pair = rtig_b + "-" + ftig_a;
									err[err_pair].links++;
									err[err_pair].gaps += d;
									ct_problem_pairs++;
									ct_problem_pairs_hash[insert_size]++;
									PET << "Pairs unsatisfied in distance within a contig pair.  rB-> <-A  WITH tig#" << tig_b << "->" << d << "<- tig#r." << tig_a << ", B=" << B_length << "  nt (start:" << B_start << ", end:" << B_end << ") A=" << A_length << " nt (start:" << A_start << ", end:" << A_end << ") CALCULATED DISTANCE APART: " << d << "<" << min_allowed << endl;
								}
							}


						}


					}
					else //Clone, paired reads located on the same contig -- could be used to investigate misassemblies
					{
					    if (verbose)
						    cout << "Pair (" << read_a << "and" << read_b << ") located on same contig " << tig_a << "(" << A_length << " nt)" << endl;
					
						int pet_size = 0;
						if (A_start > B_start && (B_start < B_end) && (A_start > A_end))
						{
							pet_size = A_start - B_start;
							track_insert[pet_size]++;
							if (pet_size >= low_iz && pet_size <= up_iz)
							{
								ct_ok_contig++;
								ct_ok_contig_hash[insert_size]++;
							}
							else
							{
								PET << "Pairs unsatisfied in distance within a contig.  Pair (" << read_a << "-" << read_b << ") on contig " << tig_a << " (" << A_length << " nt) Astart:" << A_start << " Aend:" << A_end << "Bstart:" << B_start << "Bend:" << B_end << " CALCULATED DISTANCE APART :" << pet_size << endl;
								ct_iz_issues++;
								ct_iz_issues_hash[insert_size]++;
							}
						}
						else if (B_start > A_start && (B_start > B_end) && (A_start < A_end))
						{
							pet_size = B_start - A_start;
							track_insert[pet_size]++;
							if (pet_size >= low_iz && pet_size <= up_iz)
							{
								ct_ok_contig++;
								ct_ok_contig_hash[insert_size]++;
							}
							else
							{
								PET << "Pairs unsatisfied in distance within a contig.  Pair (" << read_a << "-" << read_b << ") on contig " << tig_a << " (" << A_length << " nt) Astart:" << A_start << " Aend:" << A_end << "Bstart:" << B_start << "Bend:" << B_end << " CALCULATED DISTANCE APART :" << pet_size << endl;
								ct_iz_issues++;
								ct_iz_issues_hash[insert_size]++;
							}
						}
						else
						{
							ct_illogical++;
							ct_illogical_hash[insert_size]++;
							PET << "Pairs unsatisfied in pairing logic within a contig.  Pair (" << read_a << " -" << read_b << ") on contig" << tig_a << " (" << A_length << " nt) Astart:" << A_start << " Aend:" << A_end << " Bstart:" << B_start << "Bend:" << B_end << endl;
						}

					}

				}
				else
				{
					ct_single++;
					ct_single_hash[insert_size]++;
				}
			}//if unseen
		}//pairing read b
	}//read a

	// summary of contig pair issues
	PET << "------------- Putative issues with contig pairing - Summary  ----------------" << endl;
	vector<std::pair<string, PAIRDATA>> tmpErr;
	for (auto& ti : err)
		tmpErr.push_back(ti);
	sort(tmpErr.begin(), tmpErr.end(),
		[=](std::pair<string, PAIRDATA>& a, std::pair<string, PAIRDATA>& b) { return a.second.links > b.second.links; });
	for (int i = 0; i < tmpErr.size(); i++)
	{
		std::pair<string, PAIRDATA> te = tmpErr[i];
		string err_pair = te.first;
		double mean_iz = 0;
		if (te.second.links)
		{
			mean_iz = (double)te.second.gaps / (double)te.second.links;
		}
		PET << "Pair " << err_pair << " has " << te.second.links << " links and mean distance = " << mean_iz << endl;
	}
	PET.close();
	int satisfied = ct_ok_pairs + ct_ok_contig;
	int unsatisfied = ct_problem_pairs + ct_iz_issues + ct_illogical;
	int ct_both_reads = ct_both * 2;
	LOG << "\n===========PAIRED-END READS STATS===========" << endl;
	LOG << "At least one sequence/pair missing from contigs >= "<<contig_size_cutoff<<" bp (user-defined -z): "<<ct_single<<endl;
	LOG << "Assembled pairs: " << ct_both << "(" << ct_both_reads << " sequences)" << endl;
	LOG << "\tSatisfied in distance/logic within contigs (i.e. -> <-, distance on target: "<<ct_ok_contig<<endl;
	LOG << "\tUnsatisfied in distance within contigs (i.e. distance out-of-bounds): " << ct_iz_issues << endl;
	LOG << "\tUnsatisfied pairing logic within contigs (i.e. illogical pairing ->->, <-<- or <-->): " << ct_illogical << endl;
	LOG << "\t---" << endl;
	LOG << "\tSatisfied in distance/logic within a given contig pair (pre-scaffold): " << ct_ok_pairs << endl;
	LOG << "\tUnsatisfied in distance within a given contig pair (i.e. calculated distances out-of-bounds): " << ct_problem_pairs << endl;
	LOG << "\t---" << endl;
	LOG << "Total satisfied: " << satisfied << "\tunsatisfied: " << unsatisfied << "\n\nBreakdown by insert sizes:" << endl;

	for (auto &cIter : ct_both_hash)
	{
		int izopt = cIter.first;
		LOG << "--------Reads with " << izopt << " bp inserts--------" << endl;
		int maopt = -1 * (insert_stdev * izopt);
		int low_izopt = izopt + maopt;
		int up_izopt = izopt - maopt;
		LOG << "MIN:" << low_izopt << " MAX:" << up_izopt << " as defined by " << izopt << " * " << insert_stdev << endl;
		LOG << "At least one sequence/pair missing from contigs >=" << contig_size_cutoff << " bp (user-defined -z): " << ct_single_hash[izopt] << endl;
		LOG << "Assembled pairs: " << ct_both_hash[izopt] << endl;
		LOG << "\tSatisfied in distance/logic within contigs (i.e. -> <-, distance on target: " << ct_ok_contig_hash[izopt] << endl;
		LOG << "\tUnsatisfied in distance within contigs (i.e. distance out-of-bounds):" << ct_iz_issues_hash[izopt] << endl;
		LOG << "\tUnsatisfied pairing logic within contigs (i.e. illogical pairing ->->, <-<- or <-->): " << ct_illogical_hash[izopt] << endl;
		LOG << "\t---" << endl;
		LOG << "\tSatisfied in distance/logic within a given contig pair (pre-scaffold): " << ct_ok_pairs_hash[izopt] << endl;
		LOG << "\tUnsatisfied in distance within a given contig pair (i.e. calculated distances out-of-bounds): "<<ct_problem_pairs_hash[izopt]<<endl;

	}
	LOG << "============================================" << endl;
	ofstream CSV;
	CSV.open(distribution); 
	if (!CSV)
	{
		CSV.close();
		perror(("Can't write to " + distribution).c_str());
	}
	for (auto &tIter : track_insert)
	{
		int is = tIter.first;
		CSV << is <<","<< tIter.second << endl;

	}
	CSV.close();

}


void computeLayout(string ext, string &chain, string tig, unordered_map<string, map<string, PAIRDATA>> &pair, unordered_map<int, int> &tig_length, int &total, map<int, int>&seen_it, int orig_tig_number)
{
	string orig_tig = tig;
	int extension = 1;
	while (extension)
	{
		int tnum;
		if (regexMatch(tig,"[fr](\\d+)")&& reg_globle.matches() > 0)
			tnum = stoi(reg_globle.get_match(0));
		string tnumf = "f" + to_string(tnum);
		string tnumr = "r" + to_string(tnum);
		if (seen_it.find(tnum) == seen_it.end())
		{
			if (tnumf != orig_tig)
				seen_it[tnum]++;
			if (verbose)
				cout << "Attempt to extend tig\n";
			map<string, PAIRDATA> list = pair[tig];
			string match1 = "";
			int link1 = 0;
			int gaps1 = 0;
			string match2 = "";
			int link2 = 0;
			int gaps2 = 0;
			int cntloop = 0;

			vector<std::pair<string, PAIRDATA>> tmp;
			for (auto& ti : list)
				tmp.push_back(ti);
			sort(tmp.begin(), tmp.end(),
				[=](std::pair<string, PAIRDATA>& a, std::pair<string, PAIRDATA>& b) { return a.second.links > b.second.links; });
			//map<string, PAIRDATA> dlist;
			//list.swap(dlist);	
			for (int i = 0; i < tmp.size(); i++)
			{
				std::pair<string, PAIRDATA> t = tmp[i];
				string match = t.first;
				if (cntloop)
				{
					match2 = match;
					link2 = t.second.links;
					gaps2 = t.second.gaps;
					if (verbose)
						cout << tig<<" links second best "<<match2<<" (links:"<<link2<<" total sz:"<<gaps2<<")\n";
					break;
				}
				else
				{
					match1 = match;
					link1 = t.second.links;
					gaps1 = t.second.gaps;
					if (verbose)
						cout << tig<<" links best "<<match1<<" (links:"<<link1<<" total sz:"<<gaps1<<")\n";
				}
				cntloop++;
			}
			tmp.clear();
			tmp.shrink_to_fit();

			//###ratio
			double ratio = 0.00;
			if (link1)
				ratio = (double)link2 / (double)link1;//relative ratio of the two most abundant contig pairs
			if (regexMatch(to_string(ratio),"(\\d+\\.\\d{2})")&& reg_globle.matches() > 0)
			{
				ratio = stod(reg_globle.get_match(0));
			}
			//###mean
			int mean = 0;
			if (link1)
				mean = (int)gaps1 / (int)link1;
			int tempnum;
			if (regexMatch(match1,"[fr](\\d+)")&& reg_globle.matches() > 0)
				tempnum = stoi(reg_globle.get_match(0));
			// Assessment
			if (seen_it.find(tempnum) != seen_it.end() || link1 < min_links || ratio > max_link_ratio || tempnum == orig_tig_number)
			{
				extension = 0;
				if (verbose)
					cout << "defined seen_it->{ "<<tempnum<<" } || "<<link1 <<"<"<< min_links<<" ||"<< ratio <<">"<< max_link_ratio<<"\n L1:"<<link1<<" L2:"<<link2<<"  M1:"<<match1 <<"M2:"<<match2<<" G1:"<<gaps1<<" G2:"<<gaps2;
				break;
			}
			{
				if (verbose)
					cout << ext << "extension.  mean:" << mean << " links:" << link1 << " linkratio:" << ratio << "\n";
				if (ext == "R")
				{
					chain += "k" + to_string(link1) + "a" + to_string(ratio) + "m" + to_string(mean) + "_" + match1 + "z" + to_string(tig_length[tempnum]);
				}
				else
				{
					string temp_match = "";
					if (regexMatch(match1,"^r(\\d+)")&& reg_globle.matches() > 0)
					{
						temp_match = "f" + reg_globle.get_match(0);
					}
					else if(regexMatch(match1, "^f(\\d+)")&& reg_globle.matches() > 0)
					{
						temp_match = "r" + reg_globle.get_match(0);
					}
					chain = temp_match + "z" + to_string(tig_length[tempnum]) + "k" + to_string(link1) + "a" + to_string(ratio) + "m" + to_string(mean) + "_" + chain;
				}
				total += tig_length[tempnum];
				if (verbose)
					cout << "NEXT TIG TO LOOK AT= "<<match1<<"\n";
				tig = match1;
				extension = 1;
				if (verbose)
					cout << "Will flag "<<tnum<<" as seen  (only if "<<tnumf<<" != "<<orig_tig<<").";
				if (tnumf != orig_tig) {
					pair.erase(tnumf);
					pair.erase(tnumr);
					tig_length.erase(tnum);
				}
				else
				{
					pair.erase(tnumf);
				}
			}

		}
		else
		{
			if (verbose)
			{
				cout << "NO MORE MATCH FOR "<<tig<<" in hash: pair>>\n";
			}
			extension = 0;
			break;
		}
	}//pair is defined

}



void buildScaffolds(unordered_map<string, map<string, PAIRDATA>> &pair, unordered_map<int, int> &tig_length, int contig_size_cutoff, int verbose)
{
	map<int, int> seen_it;
	int sc_ct = 0;

	//print SC "Scaffold Number,Scaffold Size (only contig lengths considered),Scaffold Chain: e.g. _f127z7068k12a0.58m42_r3090z62k7r0.14m76_  means: contig127(+ strand=f), size 7068 (z) has 12 links (k), link ratio of 0.58 (a) and with a mean gap/overlap of 42nt (m)  with reverse (r) of contig3090 (size 62) on the right.\n";

	vector<std::pair<int, int>> tmpTl;
	for (auto& ti : tig_length)
		tmpTl.push_back(ti);
	sort(tmpTl.begin(), tmpTl.end(),
		[=](std::pair<int, int>& a, std::pair<int, int>& b) { return a.second > b.second; });
	
	for (int i = 0;i<tmpTl.size();i++)
	{
        std::pair<int, int> tlIter = tmpTl[i];
		int tig = tlIter.first;
		string ftig = "f" + to_string(tig);
		string rtig = "r" + to_string(tig);

		if (seen_it.find(tig) == seen_it.end())
		{
			//should prevent re-using a contig as seed if it's already been incorporated into a scaffold
			sc_ct++;
			string chainleft = "";

			string ori_chainright = ftig + "Z" + to_string(tlIter.second);
			string chainright = ori_chainright;
			int total = tlIter.second;

			computeLayout("R", chainright, ftig, pair, tig_length, total, seen_it, tig);
			computeLayout("L", chainleft, rtig, pair, tig_length, total, seen_it, tig);
			seen_it[tig]++;
			pair.erase(ftig);
			pair.erase(rtig);
			tig_length.erase(tig);
			string scaffold = chainleft + chainright;
			if (total >= contig_size_cutoff)
				SC << "scaffold" << sc_ct << ","<<total<<","<<scaffold << endl;
		}
		
	}
	

}


void readContigsMemory(string file, map<int, string> &fh, map<int, string> &tignames, map<int, int> &tig_length)
{
	//map<int,int> tig_length;
	//map<int,string> tignames;
	//map<int,string> fh;
	string prevhead = "";
	string seq = "";
	int cttig = 0;

	ifstream IN;
	IN.open(file);
	if (!IN)
	{
		IN.close();
		perror(("Can't read " + file).c_str());
	}
	while (!IN.eof())
	{
		string str;
		getline(IN, str);
		if (regexMatch(str,"^\\>(.*)")&& reg_globle.matches() > 0)
		{
			string head = reg_globle.get_match(0);
			if (head != prevhead && seq != "" && prevhead != "")
			{
				cttig++;
				tignames[cttig] = prevhead;
				fh[cttig] = seq;
				tig_length[cttig] = seq.length();
			}
			seq = "";
			prevhead = head;
		}
		else
		{
			transform(str.begin(), str.end(), str.begin(), ::toupper);
			seq += str;
		}
	}

	cttig++;
	tignames[cttig] = prevhead;
	fh[cttig] = seq;
	tig_length[cttig] = seq.length();

}

void buildScaffoldFasta(string dotscaffold, map<int, string>&fh, string scaffold_fasta)
{
	ifstream IN;
	IN.open(dotscaffold);
	if (!IN)
	{
		IN.close();
		perror(("Can't read " + dotscaffold).c_str());
	}

	ofstream OUT;
	OUT.open(scaffold_fasta);
	if (!OUT)
	{
		OUT.close();
		perror(("Can't write to " + scaffold_fasta).c_str());
	}

	int tot = 0;
	int ct = 0;
	int sct = 0;
	while (!IN.eof())
	{
		string str;
		getline(IN, str);
		string sc = "";
		vector<string>a = split(str, ",");
		vector<string>tig;
		if(a.size()>2)
		{
			if (regexMatch(a[2],"\\_"))
		    {
			    tig = split(a[2], "_");
		    }
		    else
		    {
			    tig.push_back(a[2]);
		    }
		}
		sct++;
		int tigsum = 0;
		if(str.size()>1)
		    OUT << ">" << str << endl;
		for (int i = 0; i < tig.size(); i++)
		{
			string t = tig[i];
			
			ct++;
			if (regexMatch(t,"([fr])(\\d+)z(\\d+)(\\S+)?","i"))
			{
				string orient = reg_globle.get_match(0);
				int tnum = stoi(reg_globle.get_match(1));
				string head = orient + to_string(tnum);
				int search = tnum;
				string other = "";
				if(reg_globle.matches() > 3)
				   other = reg_globle.get_match(3);
				tot += stoi(reg_globle.get_match(2));
				tigsum += stoi(reg_globle.get_match(2));

				string gap = "NA";
				string numlinks = "NA";
				string linksratio = "NA";

				string gapseq = "";
				if (regexMatch(other,"k(\\d+)")&& reg_globle.matches() > 0)
					numlinks = reg_globle.get_match(0);
				if (regexMatch(other, "a(\\d+.*)m")&& reg_globle.matches() > 0)
					linksratio = reg_globle.get_match(0);
				if (regexMatch(other, "m(\\-?\\d+)")&& reg_globle.matches() > 0)
					gap = reg_globle.get_match(0);
				string seq = fh[search];
				if (orient == "r")
					seq = reverseComplement(seq);
				if (gap != "NA")
				{
					if (stoi(gap) > 0)
				    {
					    for (int ii = 0; ii < stoi(gap); ii++)
					    {
						    gapseq += "N";
					    }
					}
				}	
				if (gap != "NA" && stoi(gap) <= 0)
					gapseq = "n";

				seq += gapseq;

				OUT << seq<<flush;

			}//tig regex
			
		}//each tig

		OUT << endl;

	}

	IN.close();
	OUT.close();

}


int main(int argc, char *argv[])
{
	
	//-f Herpesvirus_3.60kb.reads.fa -m 16 -o 2 -w 5 -b seedtest -c 1 -s Herpesvirus_3.60kb.seed.fa -u 1 -i 0 -j 20
	int o;
    const char *optstring = "f:m:o:v:r:p:d:l:a:z:e:g:s:t:b:c:x:n:h:u:w:j:q:y:i:"; 
    while ((o = getopt(argc, argv, optstring)) != -1) {
        switch (o) {
            case 'f':
		        opt_f = optarg;
				break;
		    case 'm':
			    opt_m = stoi(optarg);
				break;
		    case 'o':
			    opt_o = stoi(optarg);
				break;
		    case 'v':
				opt_v = stoi(optarg);
				break;
		    case 'r':
				opt_r = stoi(optarg);
				break;
		    case 'p':
				opt_p = stoi(optarg);
				break;
		    case 'd':
				opt_d = stoi(optarg);
				break;
		    case 'l':
				opt_l = stoi(optarg);
				break;
		    case 'a':
				opt_a = stoi(optarg);
				break;
		    case 'z':
				opt_z = stoi(optarg);
				break;
		    case 'e':
				opt_e = stoi(optarg);
				break;
		    case 'g':
				opt_g = optarg;
				break;
		    case 's':
				opt_s = optarg;
				break;
		    case 't':
				opt_t = stoi(optarg);
				break;
		    case 'b':
				opt_b = optarg;
				break;
		    case 'c':
				opt_c = stoi(optarg);
				break;
		    case 'x':
				opt_x = stoi(optarg);
				break;
		    case 'n':
				opt_n = stoi(optarg);
				break;
		    case 'h':
				opt_h = stoi(optarg);
				break;
		    case 'i':
				opt_i = stoi(optarg);
				break;
		    case 'w':
				opt_w = stoi(optarg);
				break;
		    case 'j':
				opt_j = stoi(optarg);
				break;
		    case 'q':
				opt_q = stoi(optarg);
				break;
		    case 'y':
				opt_y = stoi(optarg);
				break;
		    case 'u':
				opt_u = stoi(optarg);
				break;
            case '?':
                printf("error optopt: %c\n", optopt);
                printf("error opterr: %d\n", opterr);
                break;
        }
    }

	
	string file = "";
	file = opt_f;
	if (opt_m)
		min_overlap = opt_m;
	if (opt_o)
		base_overlap = opt_o;
	if (opt_r)
		min_base_ratio = opt_r;
	if (opt_t)
		max_trim = opt_t;
	if (opt_v)
		verbose = opt_v;
	if (opt_p)
		paired = opt_p;
	if (opt_l)
		min_links = opt_l;
	if (opt_a)
		max_link_ratio = opt_a;
	if (opt_z)
		contig_size_cutoff = opt_z;
	if (opt_e)
		insert_stdev = opt_e;
	if (!opt_g.empty())
		unpaired_file = opt_g;
	if (!opt_s.empty())
		seed_file = opt_s;
	if (!opt_b.empty())
		base_name = opt_b;
	if (opt_c)
		tracked = opt_c;
	if (opt_x)
		min_tig_overlap = opt_x;
	if (opt_n)
		npad_gaps = opt_n;
	if (opt_h)
		ignorehead = opt_h;
	if (opt_j)
		targetwordlen = opt_j;
	if (opt_q)
		tie_breaker = opt_q;
	if (opt_y)
		ignore_read = opt_y;
	if (opt_w)
		min_depth_of_coverage = opt_w;

	if (paired || tracked)
	{
		forcetrack = 1;
	}
	string display_unpaired_file = unpaired_file;
	string display_seed_file = seed_file;


	if (access(file.c_str(), 0) == -1)
	{
		cout << "Example:\nassembler -f paired.fa -g unpaired.fa -p 1 -m 20 -w 5 -b out" << file << endl;
		return 0;
	}
	if ((!opt_s.empty()) && (access(opt_s.c_str(), 0) == -1))//seed file specified, but does not exist
	{
		cout << "The file: " << opt_s << " you specified does not exist -- fatal" << endl;
		return 0;
	}
	else if ((!opt_s.empty()) && (access(opt_s.c_str(), 0) == 0))//seed file specified and exists
	{
		if (opt_u)
			space_restriction = opt_u;
		independent = opt_i;
		if (independent >= 1 || independent < 0)
		{
			independent = 1;
			space_restriction = 1; //restrict the space to that of the target when doing targeted de novo assembly - i mode
		}
	}

	//Naming output files
	if (base_name == "")
	{
		base_name = file + ".assbmled_m" + to_string(min_overlap) + "_o" + to_string(base_overlap) + "_r" + to_string(min_base_ratio) + "_t" + to_string(max_trim) + "_w" + to_string(min_depth_of_coverage) + "_q" + to_string(tie_breaker) + "_y" + to_string(ignore_read);


		if (paired) {
			base_name = base_name + "_e" + to_string(insert_stdev) + "_l" + to_string(min_links) + "_a" + to_string(max_link_ratio) + "_z" + to_string(contig_size_cutoff) + "_g-" + display_unpaired_file;
		}
		if (!opt_s.empty()) {
			base_name += "_s-" + display_seed_file + "_i" + to_string(independent) + "_j-" + to_string(targetwordlen) + "_u-" + to_string(space_restriction);
		}
		int pid_num = getpid();
		base_name += "_pid" + to_string(pid_num);
	}

	string contig = base_name + "_contigs.fa";
	string singlet = base_name + "_singlets.fa";
	string sshort = base_name + "_short.txt";
	string log = base_name + ".log";
	string scaffold;
	string scaffold_fasta;
	string mergedtigs;
	string issues;
	string distribution;

	if (paired)
	{
		scaffold = base_name + ".scaffolds";
		scaffold_fasta = base_name + "_scaffolds.fa";
		mergedtigs = base_name + "_mergedcontigs.fa";
		issues = base_name + "_pairing-issues.txt";
		distribution = base_name + "_pairing-distribution.csv";
	}
	string covfile;
	string rdpositionfile;
	if (tracked)
	{
		covfile = base_name + "_coverage.csv";
		rdpositionfile = base_name + ".readposition";
	}
	string pileupfile;
	if (space_restriction)
	{
		pileupfile = base_name + ".pileup";
	}


	LOG.open(log); 
	if (!LOG)
	{
		perror(("Can't write to " + log).c_str());
	}
	if (min_overlap < 16 || min_overlap > 100) 
	{
		string outofbound_message = "-m must be a number between 16-100 ...Exiting.";
		cout << outofbound_message;
		LOG << outofbound_message << endl;
		LOG.close();
		return 0;
	}

	if (base_overlap < 1) 
	{
		string outofbound_message = "-o must be set to 1 or higher ...Exiting.";
		cout << outofbound_message << endl;
		LOG << outofbound_message << endl;
		LOG.close();
		return 0;
	}

	if (min_base_ratio <= 0.5 || min_base_ratio > 1) 
	{
		string outofbound_message = "-r must be a number between 0.51 and 1.00 ...Exiting.\n";
		cout << outofbound_message << endl;
		LOG << outofbound_message << endl;
		LOG.close();
		return 0;
	}

	stringstream ss;
	ss << "\nStarting assembly reads with parameters: \n-f " << file << "\n-s " << seed_file << "\n-i " << independent << " -j " << targetwordlen << " -u " << space_restriction << " -h " << ignorehead << " -w " << min_depth_of_coverage << " -m " << min_overlap << " -o " << base_overlap << " -r " << min_base_ratio << " -t " << max_trim << " -q " << tie_breaker << " -y " << ignore_read << "\n";

	string init_message = ss.str();
	ss.clear();
	ss.str("");

	if (tracked)
	{
		ss << "-c "<<tracked<<"\n"; //Coverage: " << covfile << "\nRead position: " << rdpositionfile << "\n";
		init_message += ss.str();
		ss.clear();
		ss.str("");
		if (!independent && space_restriction)
		{
			init_message += "Pileup:" + pileupfile + "\n";
		}

	}

	if (forcetrack)
	{
		init_message += "-z " + to_string(contig_size_cutoff) + " ";
	}
	if (paired)
	{
		ss << "-p " << paired << " -e " << insert_stdev << " -l " << min_links << " -a " << max_link_ratio << "\n-g " << unpaired_file << "\n";
		init_message += ss.str();
	}
	ss.clear();
	ss.str("");
	/*init_message += "\nSinglets: " + singlet + "\nContigs: " + contig + "\n";
	if (paired)
	{
		init_message += "Scaffolds: " + scaffold_fasta + "\n";
	}
	init_message += "\nExcluded reads: " + sshort + "\nLog: " + log + "\n";
	*/
	cout << init_message;
	LOG << init_message << endl;
	assemblyruninfo = init_message + "\n";
	//---------------------------------------------
	string date = getCurrentTimeStr();
	string reading_reads_message = "\n>>>Starting reading reads: " + date + "\n";
	cout << reading_reads_message;
	LOG << reading_reads_message << endl;
	assemblyruninfo += reading_reads_message;

	unordered_map<string, string> encoded = encodeBases();
	//unordered_map<string, SEEDDATA> seed;
	unordered_map<string, int> seedsplit;

	//Sedd file in fasta formate
	if (access(opt_s.c_str(), 0) == 0)
	{
		string use_seed_message = "Using seed sequence file " + opt_s + " for this assembly.\nNote: ONLY sequences in " + opt_s + " will be used as seeds (i.e. -f " + opt_f + " and -g " + opt_g + "will NOT be used as seeds, only used for extension)";
		LOG << use_seed_message << endl;
		if (verbose)
		{
			cout << use_seed_message;
		}
		loadSeed(opt_s, targetwordlen, seed, seedsplit);

	}

	unordered_map<string, SETDATA> set;
	unordered_map<string, map<string, string>> bin;
	unordered_map<string, map<string, MATEPAIRDATA>> matepair;
	int sumall = 0;
	int ctall = 0;
	readFasta(matepair, set, bin, file, sshort, paired, encoded, seedsplit, space_restriction, targetwordlen, ignorehead, sumall, ctall);
	if (access(opt_g.c_str(), 0) == 0)
	{
		readFasta(matepair, set, bin, unpaired_file, sshort, 0, encoded, seedsplit, space_restriction, targetwordlen, ignorehead, sumall, ctall);
	}
	//unordered_map<string, int> dseedsplit;
	//seedsplit.swap(dseedsplit);

	//double fc = (double)sumall/(double)3000000000;

	string statslog = "Total reads interrogated:" + to_string(ctall) + "\tTotal base:" + to_string(sumall);
	cout << statslog << endl;
	LOG << statslog << endl;
	if (opt_s.empty())
	{
		for (auto &sIter : set)
		{
			seed[sIter.first].count = sIter.second.count;
			seed[sIter.first].names = sIter.second.names;
		}
		//seed = set;
		TRACK_COUNT = 0;
	}
	else
	{
		if (independent)
		{
			string ind_message = "-i has been set to 1, which means the target sequences are USED for recruiting reads for de novo assembly, NOT target-based reference guided assemblies\n";
			cout << ind_message;
			LOG << ind_message << endl;
			//seed = {};
			//seed = set;
			for (auto &sIter : set)
			{
				seed[sIter.first].count = sIter.second.count;
				seed[sIter.first].names = sIter.second.names;
			}
		}
		else
		{
			string ind_message = "-i has been set to 0, which means the target sequences are USED for recruiting reads for target-based reference guided assemblies (only de novo extension on the ends\n";
			cout << ind_message;
			LOG<<ind_message<<endl;
		}
	}
	string seed_number_message = "Number of unique seed sequences: " + to_string(seed.size());
	cout << seed_number_message;
	LOG << seed_number_message << endl;
	//-------------------------------------------------
	date = getCurrentTimeStr();
	string assembly_start_message = "\n>>>Assembly Initiated " + date + "\n";
	cout << assembly_start_message;
	LOG << assembly_start_message << endl;
	assemblyruninfo+=assembly_start_message + "\n"+seed_number_message+"\n";


	//-------------------------------------------------
	int sgl_count = 1;
	int tig_count = 1;
	int previous_index = 0;

	ofstream TIG;
	TIG.open(contig);
	if (!TIG)
	{
		TIG.close();
		perror(("Can't write to " + contig).c_str());
	}
	ofstream SIN;
	SIN.open(singlet); 
	if (!SIN)
	{
		SIN.close();
		perror(("Can't write to " + singlet).c_str());
	}
	
	if (paired)
	{
		SC.open(scaffold);
		if (!SC)
		{
			SC.close();
			perror(("Can't write to " + scaffold).c_str());
		}
	}
	ofstream CF;
	ofstream RP;
	if (tracked)
	{
		CF.open(covfile);
		if (!CF)
		{
			CF.close();
			perror(("Can't write to " + covfile).c_str());
		}

		RP.open(rdpositionfile);
		if (!RP)
		{
			RP.close();
			perror(("Can't write to " + rdpositionfile).c_str());
		}
	}
	ofstream PU;
	if ((!independent) && space_restriction)
	{
		PU.open(pileupfile);
		if (!PU)
		{
			PU.close();
			perror(("Can't write to " + pileupfile).c_str());
		}
	}
	unordered_map<int, int> tig_length;
	unordered_map<string, TRACKALLDATA> track_all;
	//int alternate;
	unordered_map<int, map<int, map<string,int>>> alternate;

	//try catch;
	//string status_bar = "";
	for (int i = 1; i <= 99; i++)
	{
		per[i]++;
		/*
		float ct = i / 10;
		if (ct == int(ct))
		{
			status_bar += to_string(ct);
		}
		else
		{
			status_bar += "-";
		}
		*/

	}
	//status_bar = ">>>running progress ...\n";
	//cout << status_bar << endl;
	int keys_start = seed.size();
	int low_total = 0;
	int prev_cov = 0;

	//-------------------------------------
	//#ASSEMBLY:
	vector<string>seedKeys;
	for (auto const &iter : seed)
	{
		seedKeys.push_back(iter.first);
	}
    int myi = 0;
	for (int seedi = 0; seedi < seedKeys.size(); seedi++)
	{
		string seq = seedKeys[seedi];
		unordered_map<string, TRACKDATA> track;
		string pu_seed_name = "";
		string pu_seed_ori = "";
		if (!independent && space_restriction && seed.find(seq)!=seed.end())
		{
			pu_seed_name = seed[seq].seed_name;
			pu_seed_ori = seed[seq].ori;
		}
		if (seed.find(seq) != seed.end())//#sequence read hasn't been used, is longer than 16 nt and the user-defined overlap minimum -m
		{
			//cout << myi << ":" << seq << endl;
			//myi++;
			string seed_name = "";
			if (!seed[seq].seed_name.empty())
			{
				seed_name = "|seed:" + seed[seq].seed_name;
			}

			int orig_mer = seq.length();

			if (forcetrack)
			{
				track[seq].start = 1;
				track[seq].end = orig_mer;
				track[seq].cov = seed[seq].count;
				track[seq].names = seed[seq].names;
			}

			//Delete keys ref

			string start_sequence = seq;
			int reads_needed = seed[seq].count; //tracks coverage
			int total_bases = orig_mer * reads_needed;
			deleteData(bin, set, seed, seq, encoded);//remove k-mer from hash table and prefix tree
			if (verbose)
				cout << "\n\n>>> START SEED SEQUENCE ::" << seq << "<<<\n\n";
			unordered_map<string, BARCODEDATA>  barcodeList;
			doExtension("3 prime", orig_mer, seq, set, bin, reads_needed, total_bases, min_overlap, base_overlap, min_base_ratio, verbose, track, forcetrack, tig_count, max_trim, encoded, matepair, tie_breaker, ignore_read, barcodeList);
			//end of 3' extension, beginning of 5' extension  (via 3' RC)

			string seqrc = reverseComplement(seq);
			doExtension("5 prime", orig_mer, seqrc, set, bin, reads_needed, total_bases, min_overlap, base_overlap, min_base_ratio, verbose, track, forcetrack, tig_count, max_trim, encoded, matepair, tie_breaker, ignore_read, barcodeList);
			//end of 5' extension
			//unordered_map<string, BARCODEDATA>  dbarcodeList;
			//barcodeList.swap(dbarcodeList);

			

			int leng = seqrc.length();
			string reversetig = reverseComplement(seqrc);//return to sequence, as provided 
			int trimmed_length = start_sequence.length() - 2 * max_trim;
			if (leng > trimmed_length || (reads_needed > 1 && !opt_s.empty()))
			{
				//second conditional is similar to TASR, output contig if targeted assembly mode and more reads have been recruited.
				if (low_total >= LOW_COVERAGE_TIG_IN_A_ROW)
					break;//sudden termination based on expected depth ov coverage
				double cov = (double)total_bases/(double)leng;  
				TIG << ">contig " << tig_count << "|size " << leng << "|read " << reads_needed << "|cov " << cov << seed_name << "\n" << reversetig << endl;
				//print contigs to file
				tig_length[tig_count] = leng;
				if (cov < min_depth_of_coverage)
					low_total++;
				if (forcetrack && leng >= contig_size_cutoff)
				{
					if (tracked)//only execute & report if user specifies -c
					{
						CF << ">contig " << tig_count << "|size " << leng << "|read " << reads_needed << "|cov " << cov << " " << seed_name << endl;
						RP << ">contig " << tig_count << "|size " << leng << "|read " << reads_needed << "|cov " << cov << " " << seed_name << endl;

						string tigarr = reversetig;
						map<int, PILEUPDATA> pileup;
						int gpos = 0;
						if (!independent && space_restriction)
						{
							for (int i = 0; i < tigarr.size(); i++)
							{
								string gbase;
								gbase.push_back(tigarr[i]);
								gpos++;
								pileup[gpos].tig = gbase;
							}
						}

						//initialize all positions;
						map<int, int> hashcov;
						for (int bpo = 1; bpo <= leng; bpo++)
						{
							hashcov[bpo] = 0;
						}
						vector<std::pair<string, TRACKDATA>> tmpTrack;
						for (auto& ti : track)
							tmpTrack.push_back(ti);
						sort(tmpTrack.begin(), tmpTrack.end(),
							[=](std::pair<string, TRACKDATA>& a, std::pair<string, TRACKDATA>& b) { return a.second.start < b.second.start; });
						for (int i = 0; i < tmpTrack.size(); i++)
						{
							std::pair<string, TRACKDATA> tp = tmpTrack[i];
							string rd = tp.first;
							map<string, string> rdlist = tp.second.names;
							for (auto &rItor : rdlist)
							{
								string rdoflist = rItor.first;
								if (rdoflist != "")
									RP << rdoflist << "," << tp.second.start << "," << tp.second.end << endl;
								vector<int> covposition;
								if (tp.second.start < tp.second.end)
								{
									//plus strand
									int sss = tp.second.start;
									int eee = tp.second.end;
									if (sss < 0)
										sss = 1;
									if (eee > leng)
										eee = leng;
									for (int i = sss; i <= eee; i++)
									{
										covposition.push_back(i);
									}
									string rdseq = rd;
						
									string rdqua = rItor.second;
									int rdpos = 0;
									if (!independent && space_restriction)
									{
										for (int vpos = tp.second.start; vpos <= tp.second.end; vpos++)
										{
											if (vpos >= 1 && vpos <= gpos)
											{
												if (rdoflist == pu_seed_name)
												{
													vector<string> seedarr;
													for (int pi = 0; pi < pu_seed_ori.size(); pi++)
													{
														string ts;
														ts.push_back(pu_seed_ori[pi]);
														seedarr.push_back(ts);
													}
													if (regexMatch(seedarr[rdpos],"[acgt]"))
														pileup[vpos].interest = seedarr[rdpos];


												}
												string tmpq="";
												if(rdpos<rdqua.size())
												    tmpq.push_back(rdqua[rdpos]);
												pileup[vpos].qua += tmpq;
												string trdseq="";
												if (rdpos < rdseq.size())
												    trdseq.push_back(rdseq[rdpos]);
												if (pileup[vpos].tig == trdseq)
												{
													pileup[vpos].seq += ".";
												}
												else
												{
													pileup[vpos].seq += trdseq;
												}
											}
											rdpos++;

										}
									}

								}
								else //minus strand
								{
									int sss = tp.second.end;
									int eee = tp.second.start;
									if (sss < 0)
										sss = 1;
									if (eee > leng)
										eee = leng;
									for (int i = sss; i <= eee; i++)
									{
										covposition.push_back(i);
									}
								
									string rdseq = reverseComplement(rd);
									string rstr = rdlist[rdoflist];
									reverse(rstr.begin(), rstr.end());
									string rdqua = rstr;
									int rdpos = 0;
									if (!independent && space_restriction)
									{
										for (int vpos = tp.second.end; vpos <= tp.second.start; vpos++)
										{
											if (vpos >= 1 && vpos <= gpos)
											{
												string tmpq="";
												if (rdpos < rdqua.size())
												    tmpq.push_back(rdqua[rdpos]);
												pileup[vpos].qua += tmpq;
												string trdseq="";
												if (rdpos < rdseq.size())
													trdseq.push_back(rdseq[rdpos]);
												if (pileup[vpos].tig == trdseq)
												{
													pileup[vpos].seq += ",";
												}
												else
												{
													string tr = trdseq;
													transform(tr.begin(), tr.end(), tr.begin(), ::tolower);
													pileup[vpos].seq += tr;

												}
											}
											rdpos++;
										}
									}//space restriction

								}//plus/minus


								for (int i = 0; i < covposition.size(); i++)
								{
									int pss = covposition[i];
									hashcov[pss]++;
								}


							}

							//map<string, string> drdlist;
							//rdlist.swap(drdlist);
						}
						for (auto &hIter : hashcov)
						{
							int pss = hIter.first;
							if (!hIter.second)
							{
								hIter.second = base_overlap;
							}
							CF << hIter.second << "," << flush;
							
						}
						CF << endl;
						//map<int, int> dhashcov;
						//hashcov.swap(dhashcov);
						if (!independent && space_restriction)
						{
							for (auto &pIter : pileup)
							{
								int tigpos = pIter.first;
								int depth = pIter.second.seq.length();
								string base = "";
								if (pIter.second.interest != "")
								{
									base = pIter.second.interest;
								}
								else
								{
									base = pIter.second.tig;
								}
								PU << "contig" << tig_count << "\t" << tigpos << "\t" << base << "\t" << depth << "\t" << pIter.second.seq <<" "<< pIter.second.qua << endl;
							}
							PU << endl;
						}

						//map<int, PILEUPDATA> dpileup;
						//pileup.swap(dpileup);


					}
					trackReads(track, track_all, alternate, tig_count);
					//unordered_map<string, TRACKDATA> dtrack;
					//track.swap(dtrack);
				}
				prev_cov = cov;
				tig_count++;
			}
			else
			{
				int cov = reads_needed;
				int singlet_leng = start_sequence.length();
				SIN << ">singlet " << sgl_count << "|size " << singlet_leng << "|read " << reads_needed << "|cov " << cov << " " << seed_name << "\n" << start_sequence << endl;    //print singlets to file
				sgl_count++;
			}
		}

		int keys_left = seed.size();
		double kk = (double)(keys_start-keys_left)/(double)keys_start;
		int index = (int)(kk * 100);
		if (per.find(index) != per.end())
		{
			string s = "";
			for (int i = 0; i < index - previous_index; i++)
			{
				s += ".";
			}
			//cout << s;
			per.erase(index);
		}
		previous_index = index;
		if (!keys_left)
		{
			break;
		}

	}
	//cout << ".";
	
	date = getCurrentTimeStr();
	string success = "... contig assembly completed ... " + date;
	cout << success;
	LOG << success << endl;
	assemblyruninfo += success + "\n";
	TIG.close();
	SIN.close();
	SHO.close();
	if (tracked)
	{
		CF.close();
		RP.close();
	}
	if (!independent && space_restriction) {
		PU.close();
	}

	date = getCurrentTimeStr();
	if (paired)
	{
		string sc_start_message = "\n>>>Scaffolding initiated " + date + "\n";
		cout << sc_start_message;
		LOG << sc_start_message << endl;
		assemblyruninfo += sc_start_message + "\n";
		unordered_map<string, map<string, PAIRDATA>> pair;
		pairContigs(pair, matepair, track_all, tig_length, issues, distribution, verbose);
        //unordered_map<string, map<string, MATEPAIRDATA>> dmatepair;
		//matepair.swap(dmatepair);
		//unordered_map<string, TRACKALLDATA> dtrack_all;
		//track_all.swap(dtrack_all);

		buildScaffolds(pair, tig_length, contig_size_cutoff, verbose);
		//unordered_map<int, int> dtig_length;
		//tig_length.swap(dtig_length);
		//unordered_map<string, map<string, PAIRDATA>> dpair;
		//pair.swap(dpair);

		SC.close();
		date = getCurrentTimeStr();
		string sc_end_message = "... scaffolding completed ... " + date + "\n";
		cout << sc_end_message;
		LOG << sc_end_message << endl;
		assemblyruninfo += sc_end_message + "\n";

		date = getCurrentTimeStr();
		string sc_fasta_message = "SUCCESSFULLY MAKE THE SCAFFOLD FILE: ";
		cout << sc_fasta_message;
		LOG << sc_fasta_message << endl;
		assemblyruninfo += sc_fasta_message;

		// reading contigs to make fasta file and get additional needed info
		map<int, string>tighash;
		map<int, string>tignames;
		map<int, int>tiglength;
		readContigsMemory(contig, tighash, tignames, tiglength);
		// building scaffold fasta 
		buildScaffoldFasta(scaffold, tighash, scaffold_fasta);

		string sc_fasta_done_message = scaffold_fasta+"\n";
		cout << sc_fasta_done_message;
		LOG << sc_fasta_done_message << endl;
		assemblyruninfo += sc_fasta_done_message + "\n";

	}






}

