#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <utility>
#include <algorithm>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <chrono>
// Version sin limite en memoria para obtener el score y numero de alineaciones

class NeedlemanWunschAlignment
{
private:

    std::string s;
    std::string t;

    std::vector<std::vector<bool>> up_matrix;
    std::vector<std::vector<bool>> left_matrix;
    std::vector<std::vector<bool>> diag_matrix;
    std::vector<std::vector<int>> values_matrix;

    std::vector<std::string> alignment_paths;

    int score;
    unsigned long long num_alignments;
public:
    NeedlemanWunschAlignment(std::string s, std::string t)
    {
        this->s = s;
        this->t = t;

        int rows = s.size() + 1;
        int cols = t.size() + 1;

        std::vector<int> diag_vector = {-2, -2};
        std::vector<int> prev_diag_vector = {0};
        std::vector<int> new_diag_vector;

        // Vectores para el numero de secuencias
        std::vector<unsigned long long> num_dv = {1, 1};
        std::vector<unsigned long long> num_pdv = {1};
        std::vector<unsigned long long> num_ndv;


        // Vectores para las secuencias encontradas
        std::vector<std::vector<std::string>> aligns_dv = {{"<"},{"^"}};
        std::vector<std::vector<std::string>> aligns_pdv = {{""}};
        std::vector<std::vector<std::string>> aligns_ndv;


        for (int diagonal = 2; diagonal < rows + cols - 1; diagonal++)
        {
            int start_row = std::max(0, diagonal - cols + 1);
            int start_col = std::min(diagonal, cols - 1);

            int diag_len = std::min(rows - start_row - 1, start_col) + 1;
             
            int pos_dif = 0;
            if (start_col == cols - 1 && start_row >= 2) pos_dif = 1;
            else if (start_col == cols - 1 && start_row == 1) pos_dif = 0; 
            else pos_dif = -1;

            int left_dif = 0;
            if (start_row > 0) left_dif = 1;

            int top_dif = -1;
            if (start_col == cols - 1 && start_row > 0) top_dif = 0;

            new_diag_vector.resize(diag_len);
            num_ndv = std::vector<unsigned long long>(diag_len, 0);
            aligns_ndv = std::vector<std::vector<std::string>>(diag_len);

            #pragma omp parallel for
            for (int k = 0; k < diag_len; k++)
            {
                int i = start_row +  k, j = start_col - k;
                if (i == 0)
                {
                    new_diag_vector[k] = diag_vector[k + left_dif] - 2;
                    num_ndv[k] += num_dv[k + left_dif];
                    auto tmp = aligns_dv[k + left_dif];
                    for (auto &a : tmp)
                    {
                        a += "<";
                    }
                    aligns_ndv[k].insert(aligns_ndv[k].end(), tmp.begin(), tmp.end());
                    continue;
                }
                if (j == 0)
                {
                    new_diag_vector[k] = diag_vector[k + top_dif] - 2;
                    num_ndv[k] += num_dv[k + top_dif];
                    auto tmp = aligns_dv[k + top_dif];
                    for (auto &a : tmp)
                    {
                        a += "^";
                    }
                    aligns_ndv[k].insert(aligns_ndv[k].end(), tmp.begin(), tmp.end());
                    continue;
                }

                int new_val;
                if (s[i - 1] == t[j - 1])
                {
                    new_val = std::max({diag_vector[k + top_dif] - 2, diag_vector[k + left_dif] - 2, prev_diag_vector[k + pos_dif] + 1});
                }
                else
                {
                    new_val = std::max({diag_vector[k + top_dif] - 2, diag_vector[k + left_dif] - 2, prev_diag_vector[k + pos_dif] - 1});
                }

                if (new_val == diag_vector[k + top_dif] - 2)
                {
                    num_ndv[k] += num_dv[k + top_dif];
                    auto tmp = aligns_dv[k + top_dif];
                    for (auto &a : tmp)
                    {
                        a += "^";
                    }
                    aligns_ndv[k].insert(aligns_ndv[k].end(), tmp.begin(), tmp.end());
                }
                if (new_val == diag_vector[k + left_dif] - 2)
                {
                    num_ndv[k] += num_dv[k + left_dif];
                    auto tmp = aligns_dv[k + left_dif];
                    for (auto &a : tmp)
                    {
                        a += "<";
                    }
                    aligns_ndv[k].insert(aligns_ndv[k].end(), tmp.begin(), tmp.end());
                }
                if ((new_val == prev_diag_vector[k + pos_dif] + 1 && s[i - 1] == t[j - 1]) || (new_val  == prev_diag_vector[k + pos_dif] - 1 && s[i - 1] != t[j - 1]))
                {
                    num_ndv[k] += num_pdv[k + pos_dif];
                    auto tmp = aligns_pdv[k + pos_dif];
                    for (auto &a : tmp)
                    {
                        a += "\\";
                    }
                    aligns_ndv[k].insert(aligns_ndv[k].end(), tmp.begin(), tmp.end());
                }

                new_diag_vector[k] = new_val;
            }
            prev_diag_vector = diag_vector;
            diag_vector = new_diag_vector;
            new_diag_vector.clear();

            num_pdv = num_dv;
            num_dv = num_ndv;
            num_ndv.clear();

            aligns_pdv = aligns_dv;
            aligns_dv = aligns_ndv;
            aligns_ndv.clear();
        }

        this->score = diag_vector[0];
        this->num_alignments = num_dv[0];
        this->alignment_paths = aligns_dv[0];
    }

    int get_score()
    {
        // return this->values_matrix[s.size()][t.size()];
        return this->score;
    }

    std::vector<std::pair<std::string, std::string>> get_allignments()
    {
        return this->get_allignments(this->s.size(), this->t.size());
    }
    std::vector<std::pair<std::string, std::string>> get_allignments_parallel()
    {
        std::vector<std::pair<std::string, std::string>> res;

        for (int i = 0; i < alignment_paths.size(); i++)
        {
            int x = s.size(), y = t.size();
            std::string s_suffix;
            std::string t_suffix;
            int j = alignment_paths[i].size() - 1;
            while (x > 0 || y > 0)
            {
                if (alignment_paths[i][j] == '<')
                {
                    s_suffix.insert(0, 1, '_');
                    t_suffix.insert(0, 1, this->t[y - 1]);
                    y--;
                }
                else if (alignment_paths[i][j] == '^')
                {
                    s_suffix.insert(0, 1, this->s[x - 1]);
                    t_suffix.insert(0, 1, '_');
                    x--;
                }
                else 
                {
                    s_suffix.insert(0, 1, this->s[x - 1]);
                    t_suffix.insert(0, 1, this->t[y - 1]);
                    x--;
                    y--;
                }
                j--;
            }
            res.push_back(std::make_pair(s_suffix, t_suffix));
        }
        return res;
    }

    unsigned long long optimal_allignments_num()
    {
        return this->num_alignments;
    }

    std::pair<std::string, std::string> get_one_alignment()
    {
        int i = this->s.size();
        int j = this->t.size();

        std::string alignment_s = "";
        std::string alignment_t = "";

        while (i > 0 || j > 0)
        {
            if (i == 0)
            {
                while (j > 0)
                {
                    alignment_s.insert(0, 1, '_');
                    alignment_t.insert(0, 1, this->t[j - 1]);
                    j--;
                }
                break;
            }
            else if (j == 0)
            {
                while (i > 0)
                {
                    alignment_s.insert(0, 1, this->s[i - 1]);
                    alignment_t.insert(0, 1, '_');
                    i--;
                }
                break;
            }
            else if (diag_matrix[i][j])
            { 
                alignment_s.insert(0, 1, this->s[i - 1]);
                alignment_t.insert(0, 1, this->t[j - 1]);
                i--;
                j--;
            }
            else if (up_matrix[i][j])
            {
                alignment_s.insert(0, 1, this->s[i - 1]);
                alignment_t.insert(0, 1, '_');
                i--;
            }
            else
            {
                alignment_s.insert(0, 1, '_');
                alignment_t.insert(0, 1, this->t[j - 1]);
                j--;
            }
        }

        return std::make_pair(alignment_s, alignment_t);
    }

    std::vector<std::pair<std::string, std::string>> get_allignments(int i, int j)
    {
        std::vector<std::pair<std::string, std::string>> allignments;
        std::string s_suffix;
        std::string t_suffix;

        while (i > 0 || j > 0)
        {
            if (i == 0)
            {
                while (j > 0)
                {
                    s_suffix.insert(0, 1, '_');
                    t_suffix.insert(0, 1, this->t[j - 1]);
                    j--;
                }
                break;
            }
            else if (j == 0)
            {
                while (i > 0)
                {
                    s_suffix.insert(0, 1, this->s[i - 1]);
                    t_suffix.insert(0, 1, '_');
                    i--;
                }
                break;
            }
            else if (diag_matrix[i][j] && !up_matrix[i][j] && !left_matrix[i][j])
            { 
                s_suffix.insert(0, 1, this->s[i - 1]);
                t_suffix.insert(0, 1, this->t[j - 1]);
                i--;
                j--;
            }
            else if (!diag_matrix[i][j] && up_matrix[i][j] && !left_matrix[i][i])
            {
                s_suffix.insert(0, 1, this->s[i - 1]);
                t_suffix.insert(0, 1, '_');
                i--;
            }
            else if (!diag_matrix[i][j] && !up_matrix[i][j] && left_matrix[i][i])
            {
                s_suffix.insert(0, 1, '_');
                t_suffix.insert(0, 1, this->t[j - 1]);
                j--;
            }
            else
            {
                if (diag_matrix[i][j])
                {
                    std::vector<std::pair<std::string, std::string>> new_allign = this->get_allignments(i - 1, j - 1);
                    for (int k = 0; k < new_allign.size(); k++)
                    {
                        allignments.push_back(std::make_pair(std::string(new_allign[k].first + this->s[i - 1] + s_suffix), std::string(new_allign[k].second + this->t[j - 1] + t_suffix)));
                    }
                }
                if (up_matrix[i][j])
                {
                    std::vector<std::pair<std::string, std::string>> new_allign = this->get_allignments(i - 1, j);
                    for (int k = 0; k < new_allign.size(); k++)
                    {
                        allignments.push_back(std::make_pair(std::string(new_allign[k].first + this->s[i - 1] + s_suffix), std::string(new_allign[k].second + '_' + t_suffix)));
                    }
                }
                if (left_matrix[i][j])
                {
                    std::vector<std::pair<std::string, std::string>> new_allign = this->get_allignments(i, j - 1);
                    for (int k = 0; k < new_allign.size(); k++)
                    {
                        allignments.push_back(std::make_pair(std::string(new_allign[k].first + '_' + s_suffix), std::string(new_allign[k].second + this->t[j - 1] + t_suffix)));
                    }
                }
                break;
            }
        }

        // Si solo encontro un camino
        if (allignments.size() == 0)
        {
            allignments.push_back(std::make_pair(s_suffix, t_suffix));
        }
        return allignments;
    }


};

std::string get_random_dna_sequence(int length)
{
    std::string alphabet = "ACGT";
    std::string random_seq = "";
    for (int i = 0; i < length; i++)
    {
        int p = rand() % 4;
        random_seq.push_back(alphabet[p]);
    }
    return random_seq;
}

std::pair<std::string, std::string> get_sequences_from_file(std::string filename)
{
    std::ifstream file(filename.c_str());
    std::string s, t;
    file >> s;
    file >> t;

    file.close();

    return std::make_pair(s, t); 
}



int main()
{
    srand(time(NULL));

    std::pair<std::string, std::string> dna = get_sequences_from_file("dna.txt");

    std::string s = get_random_dna_sequence(10);
    std::string t = get_random_dna_sequence(10);

    std::ofstream file("result_opt.txt");
    const auto start{std::chrono::steady_clock::now()};
    NeedlemanWunschAlignment alineacion(s, t);
    const auto end{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_seconds{end - start};
    std::cout << elapsed_seconds.count() << '\n';

    std::vector<std::pair<std::string, std::string>> found_alignments = alineacion.get_allignments_parallel();


    file << "Resultados para las cadenas:\n";
    file << "S: " << s << '\n';
    file << "T: " << t << "\n\n";
    file << "Score: " << alineacion.get_score() << '\n';

    file << "Cantidad de alineamientos producidos: " << alineacion.optimal_allignments_num() << '\n';
    file << "Alineamientos generados:\n";


    for (int i = 0; i < found_alignments.size(); i++)
    {
        file << found_alignments[i].first << '\n' << found_alignments[i].second << "\n\n";
    }

    file.close();

    return 0;
}