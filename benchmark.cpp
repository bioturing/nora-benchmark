#include <string>
#include <map>
#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <set>

#include <malloc.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#define LOG_OFFSET		0.01
#define EPSILON			1e-300
#define DIGIT			100	/* round up to 2 digit after decimal */

void get_rank(std::vector<double> *rank,
	      const std::vector<double> &tpm)
{
	std::vector<std::pair<double, int> > temp;
	temp.resize(tpm.size());
	int i, j, n = (int)temp.size();
	
	for (i = 0; i < n; i++)
		temp[i] = std::make_pair(tpm[i], i);
	sort(temp.begin(), temp.end());
	(*rank).resize(temp.size());

	for (i = 0; i < n; ) {
		double val, r;
		
		val = temp[i].first;
		j = i + 1;
		while (j < n && val == temp[j].first)
			j++;

		r = (i + j - 1) / 2.0;
		for ( ; i < j; i++)
			(*rank)[temp[i].second] = r;
	}
}

void parse_truth(char *file_path,
		 std::vector<double> *truth,
		 std::map<std::string, int> *truth_map)
{ 
	std::ifstream fi(file_path);
	if (!fi.is_open()) {
		fprintf(stderr, "ERROR: Could not open truth file!\n");
		exit(EXIT_FAILURE);
	}

	std::string name;
	double tpm;
	int id = 0;

	while (fi >> name >> tpm) {
		tpm = round(tpm * DIGIT) / DIGIT;
		(*truth).push_back(tpm);
		if ((*truth_map).count(name)) {
			fprintf(stderr, "ERROR: Duplicated transcripts %s in truth file!\n",
				name.c_str());
			exit(EXIT_FAILURE);
		}
		(*truth_map)[name] = id++;
	}
}

void parse_input(char *file_path,
		std::vector<double> *input,
		std::map<std::string, int> &truth_map)
{
	int id, num;
	std::string name;
	std::set<std::string> p;
	double tpm;

	std::ifstream fi(file_path);
	if (!fi.is_open()) {
		fprintf(stderr, "ERROR: Could not open truth file!\n");
		exit(EXIT_FAILURE);
	}

	(*input).resize(truth_map.size());
	num = 0;

	while (fi >> name >> tpm) {
		if (!truth_map.count(name)) {
			fprintf(stderr, "ERROR: Could not find corresponding transcript %s in truth file!\n",
				name.c_str());
			exit(EXIT_FAILURE);
		}
		if (p.count(name)) {
			fprintf(stderr, "ERROR: Duplicated transcripts %s in input file!\n",
				name.c_str());
			exit(EXIT_FAILURE);
		}
		p.insert(name);
		id = truth_map[name];
		tpm = round(tpm * DIGIT) / DIGIT;
		(*input)[id] = tpm;
		++num;
	}

	if (num != (int)truth_map.size()) {
		fprintf(stderr, "ERROR: Number of transcripts in two files are not equal!\n");
		exit(EXIT_FAILURE);
	}

	p.clear();
}

double mae(const std::vector<double> &truth, const std::vector<double> &input, int n)
{
	double mean = 0;
	int i, non_zero = 0;

	for (i = 0; i < n; i++) {
		double diff = abs(asinh(truth[i]) - asinh(input[i]));
		non_zero += (truth[i] > EPSILON || input[i] > EPSILON);
		mean += diff;
	}

	assert(non_zero != 0);
	mean /= non_zero;
	return mean;
}

double pearson(const std::vector<double> &truth, const std::vector<double> &input, int n)
{
	double mean1 = 0, mean2 = 0;
	double var1 = 0, var2 = 0, num = 0;
	int i = 0;

	for (i = 0; i < n; i++) {
		mean1 += truth[i];
		mean2 += input[i];
	}
	
	mean1 /= n;
	mean2 /= n;

	for (i = 0; i < n; i++) {
		num += (truth[i] - mean1) * (input[i] - mean2);
		var1 += (truth[i] - mean1) * (truth[i] - mean1);
		var2 += (input[i] - mean2) * (input[i] - mean2);
	}

	return num / sqrt(var1 * var2);
}

double log_pearson(const std::vector<double> &truth, const std::vector<double> &input, int n)
{
	double mean1 = 0, mean2 = 0;
	double var1 = 0, var2 = 0, num = 0;
	int i = 0;

	for (i = 0; i < n; i++) {
		mean1 += log(truth[i] + LOG_OFFSET);
		mean2 += log(input[i] + LOG_OFFSET);
	}
	
	mean1 /= n;
	mean2 /= n;

	for (i = 0; i < n; i++) {
		double x = log(truth[i] + LOG_OFFSET) - mean1;
		double y = log(input[i] + LOG_OFFSET) - mean2;
		num += x * y;
		var1 += x * x;
		var2 += y * y;
	}

	return num / sqrt(var1 * var2);
}

void scoring(const std::vector<double> &t_rank, const std::vector<double> &i_rank,
	     const std::vector<double> &truth, const std::vector<double> &input)
{
	int max_fp, max_fn, n_fn, n_fp, i, n = truth.size();
	double tpm_fn = 0, tpm_fp = 0;

	fprintf(stdout, "Spearman:\t\t\t%.6lf\n", pearson(t_rank, i_rank, n));
	fprintf(stdout, "Pearson:\t\t\t%.6lf\n", pearson(truth, input, n));
	fprintf(stdout, "Log-pearson:\t\t\t%.6lf\n", log_pearson(truth, input, n));
	fprintf(stdout, "MAE(asinh):\t\t\t%.6lf\n", mae(truth, input, n));	

	max_fp = -1, max_fn = -1, n_fn = 0, n_fp = 0;
	for (i = 0; i < n; i++) {
		/* false negative */
		if (truth[i] > EPSILON && input[i] <= EPSILON) {	
			if (max_fn < 0 || truth[i] > truth[max_fn])
				max_fn = i;
			++n_fn;
			tpm_fn += truth[i];
		/* false positive */
		} else if (truth[i] <= EPSILON && input[i] > EPSILON) {	
			if (max_fp < 0 || input[i] > input[max_fp])
				max_fp = i;
			++n_fp;
			tpm_fp += input[i];
		}
	}

	fprintf(stdout, "False negative:\t\t\t%d (%.2f%%)\n", n_fn, 100.0 * n_fn / n);
	fprintf(stdout, "False positive:\t\t\t%d (%.2f%%)\n", n_fp, 100.0 * n_fp / n);
	fprintf(stdout, "Max false negative:\t\t%lg\n", truth[max_fn]);
	fprintf(stdout, "Max false positive:\t\t%lg\n", input[max_fp]);
	fprintf(stdout, "Total false negative tpm:\t%lg\n", tpm_fn);
	fprintf(stdout, "Total false positive tpm:\t%lg\n", tpm_fp);
}

void print_usage()
{
	fprintf(stderr, "Developed by BioTuring (www.bioturing.com), this tool calculates the correlation between two\n");
	fprintf(stderr, "transcript expression files. For more details, please check the README on Github:\n");
	fprintf(stderr, "https://github.com/bioturing/nora-benchmark\n"); 
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage: ./benchmark TRUTH_FILE INPUT_FILE\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Truth file and input must be a tsv (tab-separated) file — with no header. Containing two\n");
	fprintf(stderr, "columns mapping of each transcript present in the reference to the corresponding tpm values\n");
	fprintf(stderr, "(the first column is a transcript and the second is the corresponding tpm value).\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "ENST00000000233    28.410000\n");
	fprintf(stderr, "ENST00000000412    0.000000\n");
	fprintf(stderr, "ENST00000000442    0.000000\n");
	fprintf(stderr, "...\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "You can use this command to convert Nora output file to this format:\n");
	fprintf(stderr, "   tail -n +2 nora_out.tsv | awk '{printf \"%%s\\t%%s\\n\", $1, $8}' > input.tsv\n");
	fprintf(stderr, "\n");
}

int main(int argc, char **argv)
{
	if (argc != 3) {
		print_usage();                      
		return 1;
	}

	std::vector<double> truth, input, i_rank, t_rank;
	std::map<std::string, int> truth_map;

	parse_truth(argv[1], &truth, &truth_map);
	parse_input(argv[2], &input, truth_map);
	get_rank(&t_rank, truth);
	get_rank(&i_rank, input);
	scoring(t_rank, i_rank, truth, input);

        return 0;
}
