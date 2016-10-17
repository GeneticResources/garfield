// Copyright (C) 2014 Genome Research Ltd / EMBL - European Bioinformatics Institute
//
// Author : Valentina Iotchkova <vi1@sanger.ac.uk>
// Author : Matthias Geihs <mg18@sanger.ac.uk>
// Author : Sandro Morganella <sm22@sanger.ac.uk>
//
// This file is part of the package garfield - GWAS analysis of regulatory or functional information
// enrichment with LD correction.
//
// GARFIELD is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 3
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, see <http://www.gnu.org/licenses/>.


#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <R.h>

using namespace std;

int n;
vector<int> vnp, vna;
vector<vector<int> > vvnap;
ofstream myfile;
// timing

clock_t tic() {
	return clock();
}

clock_t toc(clock_t t, const string& a = string()) {
	Rprintf("\t\t\t");
	//cerr << "Time elapsed: " << 1.0 * (clock() - t) / CLOCKS_PER_SEC << "s";
	Rprintf("Time elapsed: %f s", 1.0 * (clock() - t) / CLOCKS_PER_SEC);
	//string astr = a;
	if (!a.empty()) Rprintf(" (%s)",a.c_str());
	Rprintf("\n");
	return clock();
}


// bin container
class bin {
private:
	vector<double> p;
	vector< vector<bool> > a;
public:
	void put(double p, vector<bool> a) {
		this->p.push_back(p);
		this->a.push_back(a);
	}
	double getp(long i) const {
		return this->p[i];
	}
	bool geta(long i, long ia) const {
		return this->a[i][ia];
	}
	long size() const {
		return this->p.size();
	}
};

// compute n adjusted quantiles of v and store in q
void quantiles(vector<double>& q, const vector<double>& v, int n) {
	vector<double> v2 = v;
	sort(v2.begin(), v2.end());
	long vstart = 0;
	q.clear();
	for (int i=n; i>0; i--) {
		if (vstart+1 > ( (long) v2.size())) {
			Rprintf("Warning: Created %d quantiles (%d requested).\n", n-i, n);
			break;
		}
		long j = vstart + (v2.size() - vstart - 1) / i;
		q.push_back(v2[j]);
		for (; j<( (long) v2.size()); j++) {
			if (v2[j] > q.back()) break;
		}
		vstart = j;
	}
}

// get quantile index of x having quantiles q
long qindex(const vector<double>& q, double x) {
	for (long i=0; i<( (long) q.size()); i++) {
		if (x <= q[i]) return i;
	}
	throw runtime_error("Error: Could not find quantile, value out of range!");
}

// count annotated snps in bin for annotation ia
long countna(const bin& b, long ia) {
	long na = 0;
	for (long i=0; i<b.size(); i++) {
		if (b.geta(i,ia)) na++;
	}
	return na;
}

// count annotated snps with pval<=pt in all bins for annotation ia
long countnap(const vector<bin>& vbin, double pt, long ia) {
	long nap = 0;
	for(vector<bin>::const_iterator b = vbin.begin(); b != vbin.end(); b++) {
		for (long i=0; i<b->size(); i++) {
			if (b->getp(i)<=pt && b->geta(i,ia)) nap++;
		}
	}
	return nap;
}

// count snps with pval<=pt in bin
long countnp(const bin& b, double pt) {
	long np = 0;
	for (long i=0; i<b.size(); i++) {
		if (b.getp(i) <= pt) np++;
	}
	return np;
}

// vd = vectors of probability distribution for all bins, nperm = number of permutations, nbins = number of bins, nap = observed number of annotated variants at threshold, optim_mode = true if optim_mode algorith to be used, thresh = threshold for significance testing if optim_mode algorithm is used, minit = min number of iterations if optim_mode algorithm is used 
// simulate fold enrichment values and calculate empirical p-value distribution
double ep(vector< vector<double> > vd, long nperm, long nbins, int nap, int optim_mode, double thresh, int minit){

	double randunif, epval=0.0;
	int randnap[nperm], rdd, ind=0;
	
	//create an array of nperm by nbins random numbers 
	for (int i=0; i<nperm; i++){
		randnap[i]=0;
		for (int j=0; j<nbins; j++){

			randunif = (rand()+0.0)/(RAND_MAX+0.0);
		
			// here take vd[j] and create nap random numbers
			rdd = 0;
			if (vd[j].size()>1){
				for (long ii=1; ii<( (long) vd[j].size()); ii++){ // can sort and save some cpu time...
					if (randunif>=vd[j][ii-1]){
						rdd = rdd+1;
					} else {
						break;
					}	
				}
			}
			randnap[i] += rdd;
		}
		if (randnap[i]>=nap) epval+=1.0;

		// option for optim_mode iterative algorithm
		if (optim_mode == 1){
			if (epval>= thresh*nperm && i>=minit){
				ind=i+1;
				epval /= (i+1);
				break;
			}
		}
	}
	if (ind == 0){
		epval /= nperm;
	}
	return(epval);
}

// n = population size, np = number of successes, na = number of trials
//discrete_distribution<int> generate_hg_distribution(long n, long np, long na) {
vector<double> generate_hg_distribution(long n, long np, long na) {
	if (na<=0 || n<na) return vector<double>();//if (na<=0 || n<na) return discrete_distribution<int>();

	long minnanp = min(np, na);
	vector<double> v(minnanp+1), vv(minnanp+1);

	{
		double p = 0.0;
		for (long i=0; i<na; i++) {
			p += log(n-np-i) - log(n-i);
		}
		v[0] = exp(p);
		vv[0] = v[0];
	}
	
	for (long i=1; i<=minnanp; i++) {
		v[i] = v[i-1] * (np - i + 1) * (na - i +1) / (i * ((n - np) - na + i));
		vv[i] = vv[i-1] + v[i];
	}

	//return discrete_distribution<int>(v.begin(), v.end());
	return vv;//discrete_dist(v);
}



// compute fold enrichment and empirical pvalues
void batchfeep(vector< vector<double> >& vvfe, vector< vector<double> >& vvep,
	const vector<bin>& vbin, const vector<double>& vpt,
	int nannot, int nperm, bool show_progress, int optim_mode, double thresh, int minit, vector<bool> test) {
	clock_t t = tic();

	// precomputations

	// precompute na for all bins and annotations
	vna.resize(nannot);
	vector< vector<long> > vvna(nannot, vector<long>(vbin.size()));
	for (int i=0; i<nannot; i++) {
		long nasum = 0;
		for(long j=0; j<( (long) vbin.size()); j++) {
			long na = countna(vbin[j], i);
			vvna[i][j] = na;
			nasum += na;
		}
		vna[i] = nasum;
	}

	// precompute nap for all pval thresholds and annotations
	vvnap.resize(vpt.size());

	for (long i=0; i<( (long) vpt.size()); i++) {
		vvnap[i].resize(nannot);
		for (long j=0; j<nannot; j++) {
			vvnap[i][j] = countnap(vbin, vpt[i], j);
		}
	}

	// store snp count
	n = 0;
	for(long i=0; i<( (long) vbin.size()); i++) {
		n += vbin[i].size();
	}

	// Valentina added global variables to be output
	//vna_global = vna;
	//vvnap_global = vvnap;

	// precompute np for all bins and pval thresholds
	vnp.resize(vpt.size());
	vector< vector<long> > vvnp(vpt.size(), vector<long>(vbin.size()));
	for (int i=0; i<( (long) vpt.size()); i++) {
		long npsum = 0;
		for(long j=0; j<( (long) vbin.size()); j++) {
			long np = countnp(vbin[j], vpt[i]);
			vvnp[i][j] = np;
			npsum += np;
		}
		vnp[i] = npsum;
	}

	// precompute hypergeometric distributions for pval thresholds, annotations and bins
	//vector< vector< vector< discrete_distribution<int> > > > vvvd(vpt.size(), vector< vector< discrete_distribution<int> > >(nannot, vector< discrete_distribution<int> >(vbin.size())));
	vector< vector< vector <vector<double> > > > vvvd(vpt.size(), vector< vector< vector<double> > >(nannot, vector< vector<double> >(vbin.size())));
	for (long i=0; i<( (long) vpt.size()); i++) {
		for (long j=0; j<nannot; j++) {
			for (long k=0; k<( (long) vbin.size()); k++) {
				vvvd[i][j][k] = generate_hg_distribution(vbin[k].size(), vvnp[i][k], vvna[j][k]);
			}
		}
	}

	t = toc(t, "precalc");

	
	
	// compute fold enrichment for each pval threshold and annotation
	vvfe = vector< vector<double> >(vpt.size(), vector<double>(nannot));
	for (long i=0; i<( (long) vpt.size()); i++) {
		for (long j=0; j<nannot; j++) {
			if(vnp[i]==0 || vna[j]==0){
				vvfe[i][j]=-1;
			}else{
				vvfe[i][j] = (1.0 * vvnap[i][j] / vna[j]) / (1.0 * vnp[i] / n);
			}
		}
	}
	// compute empirical pvalue for each pval threshold and annotation simulating nperm permutations
	srand (time(NULL));	
	vvep = vector< vector<double> >(vpt.size(), vector<double>(nannot));
	
	clock_t t0 = tic();
	for (long j=0; j<nannot; j++) {
		for (long i=0; i<( (long) vpt.size()); i++) {
			if (test[i]) {
				vvep[i][j] = ep(vvvd[i][j], nperm, vbin.size(), vvnap[i][j], optim_mode, thresh, minit);
			} else {
				vvep[i][j] = -1;
			}
		}
		if (show_progress) {
		float perms_per_sec = 1.0 * (j+1) / (1.0 * (clock() - t0) / CLOCKS_PER_SEC);
		Rprintf("\r%s/%s\t\t%s\%\t\t%s annotations/s\t\t Time remaining: %s s",j+1,nannot, 100.0 * (j+1) / nannot, nannot,1.0 * (nannot - (j+1)) / perms_per_sec   );
		//cerr << '\r' << j+1 << '/' << nannot << "    " << 100.0 * (j+1) / nannot << '%'
		//	 << "    " << nannot << " annotations/s"
		//	 << "    " << "Time remaining: " << 1.0 * (nannot - (j+1)) / perms_per_sec << 's';
		}
	}
	if (show_progress) Rprintf("\n");
	t = toc(t, "calc");
}

extern "C" int garfield_perm(char **input, char **link, char **out_file, char **p_thresh, char **pt_thresh, char **npermut_par, char **nannot_par, char **nqmaf_par, char **nqntag_par, char **nqtssd_par, char **optim_mode_par, char **minit_par, char **thresh) {
	
	int npermut = atoi((string(*npermut_par)).c_str());
	int nannot = atoi((string(*nannot_par)).c_str());
	int nqmaf = atoi((string(*nqmaf_par)).c_str());
	int nqntag = atoi((string(*nqntag_par)).c_str());
	int nqtssd = atoi((string(*nqtssd_par)).c_str());
	int optim_mode = atoi((string(*optim_mode_par)).c_str());
	int minit = atoi((string(*minit_par)).c_str());
	double thresh0 = atof((string(*thresh)).c_str());
	/*Rprintf("qui in> %s\n", input.c_str());
	Rprintf("qui link> %s\n", link.c_str());
	Rprintf("qui out> %s\n", out_file.c_str());
	Rprintf("qui p_thr> %s\n", p_thresh.c_str());
	Rprintf("qui pt_thr> %s\n", pt_thresh.c_str());
	Rprintf("qui nannot>%d npermut>%d maf>%d tag>%d tssd>%d optim>%d mint>%d thresh>%f", npermut, nannot, nqmaf, nqntag, nqtssd, optim_mode, minit, thresh);
	*/
	clock_t tstart = tic();
	clock_t t = tic();
	
   	vector<double> pthresh, pthreshtest;
   	
   	stringstream ss1((string(*p_thresh)).c_str());
	string x1;
	while (getline(ss1, x1, ',')) {
	pthresh.push_back(atof(x1.c_str()));
	}
	
	stringstream ss2((string(*pt_thresh)).c_str());
	string x2;
	while (getline(ss2, x2, ',')) {
		pthreshtest.push_back(atof(x2.c_str()));
	}
				
	if (npermut <= 0) {
		Rprintf("\t\t[ERROR] Number of permutations must be greater than 0!\n");
		return(1);
	}
	
	if (nannot <= 0) {
		Rprintf("%s",nannot);
		Rprintf("\t\t[ERROR] Number of annotations must be greater than 0!\n");
        return(1);
	}
	
	if (pthresh.size() == 0) {
		Rprintf("\t\t[ERROR] Please provide at least one pval threshold for fold enrichment calculation!\n");
        return(1);
	}
	
	if (pthreshtest.size() == 0) {
		Rprintf("\t\t[ERROR] Please provide at least one pval threshold for significance testing!\n");
        return(1);
	}

	if (nqntag == 0 || nqmaf == 0 || nqtssd == 0) {
		Rprintf("\t\t[ERROR] Number of quantiles cannot be 0!\n");
        return(1);
	}
	
	if (string(*input).compare("NULL")==0) {
		Rprintf("\t\t[ERROR] No input filename given!\n");
        return(1);
    	}
	int ind=0;
	vector<bool> test(pthresh.size(), false);
	for (int i=0; i<( (int) test.size()); i++){
		for (int j=0; j<( (int) pthreshtest.size()); j++){
			if (pthresh[i] == pthreshtest[j]) {
				test[i] = true;
				ind += 1;
			}
		}
	}
	if (ind != ( (int) pthreshtest.size())){
		Rprintf("Not all pvalue thresholds to be tested belong to those with fold enrichment calculation!\n");
		return(1);
	}	
	// check that pthresh contains all values of pthreshtest
	
	
	ofstream out((string(*out_file)).c_str());
	if (!out.is_open()) {
		Rprintf("\t\t[ERROR]  Could not open output file!\n");
	       	return(1);
	}
	// compute binbounds and create bins
	
	
	vector<double> qntag, qmaf, qtssd;
	vector<bin> bins;
	string line;
	
	
	string fn((string(*input)).c_str());
	ifstream f(fn.c_str());
	if (!f.is_open()) {
		Rprintf("Could not open input file!\n");
	        return(1);
	}
	vector<double> vntag, vmaf, vtssd;
	vector<double> vp;
	vector< vector<bool> > vva;
	  
	while (getline(f, line)) {
	    	// snpid pval ntag maf tssd a1a2...
	    	istringstream iss(line);
	    	string snpid; double pval; long ntag; double maf; long tssd; string a;

	    	iss >> snpid >> pval >> ntag >> maf >> tssd >> a;
	    	if (iss.fail() || ((int) a.length()) < nannot) throw runtime_error("Input file has bad format!");

	    	vp.push_back(pval);
	    	vntag.push_back(ntag);
			vmaf.push_back(maf);
			vtssd.push_back(tssd);

			vector<bool> va(nannot, false);
    		for (int i=0; i<( (int) va.size()); i++) {
    			va[i] = (a[i]=='1');
    		}
    		vva.push_back(va);
	}
	f.close();

	// read link file
	string ln((string(*link)).c_str());
	ifstream l(ln.c_str());
	if (!l.is_open()) {
		Rprintf("Could not open link file!\n");
	        return(1);
	}
	vector< vector<string> > lnk;
	while (getline(l, line)) {
	    	// snpid pval ntag maf tssd a1a2...
	    	istringstream iss(line);
	    	string cl1, cl2, cl3, cl4, cl5, cl6;

		vector<string> tmp(6);
    		for (int i=0; i<( (int) tmp.size()); i++) {
    			iss >> tmp[i];
    		}
    		lnk.push_back(tmp);
	}
	l.close();



	t = toc(t, "read file");
		
	quantiles(qntag, vntag, nqntag); nqntag = qntag.size();
	quantiles(qmaf, vmaf, nqmaf); nqmaf = qmaf.size();
	quantiles(qtssd, vtssd, nqtssd); nqtssd = qtssd.size();

	bins = vector<bin>(nqntag * nqmaf * nqtssd);

	for (unsigned long i=0; i<vp.size(); i++)
		bins[qindex(qmaf, vmaf[i]) +
		qindex(qntag, vntag[i]) * nqmaf +
		qindex(qtssd, vtssd[i]) * (nqmaf * nqntag)].put(vp[i], vva[i]);

	t = toc(t, "create bins");
	

	// calculate fold enrichment and empirical pvalue
	vector< vector<double> > vvfe, vvep;
	batchfeep(vvfe, vvep, bins, pthresh, nannot, npermut, FALSE, optim_mode, thresh0, minit, test);
	t = toc(t, "calculate fold enrichment and empirical pvalue");
	out<<"Index PThresh FE EmpPval NAnnotThesh NAnnot NThresh N "<<lnk[0][1]<<" "<<lnk[0][2]<<" "<<lnk[0][3]<<" "<<lnk[0][4]<<" "<<lnk[0][5]<<endl;	
	for (int i=0; i<nannot; i++) {
		for (long j=0; j<( (long) pthresh.size()); j++) {
			out<<i<<" "<<pthresh[j]<<" "<<vvfe[j][i]<<" "<<vvep[j][i]<<" "<<vvnap[j][i]<<" "<<vna[i]<<" "<<vnp[j]<<" "<<n<<" "<<lnk[i+1][1]<<" "<<lnk[i+1][2]<<" "<<lnk[i+1][3]<<" "<<lnk[i+1][4]<<" "<<lnk[i+1][5]<<endl;
		}
	}

	toc(tstart, "total");
	out.close();
	return(0);
}
