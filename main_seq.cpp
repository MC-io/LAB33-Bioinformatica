#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <utility>
#include <algorithm>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <chrono>

class NeedlemanWunschAlignment
{
public:

    std::string s;
    std::string t;

    std::vector<std::vector<bool>> up_matrix;
    std::vector<std::vector<bool>> left_matrix;
    std::vector<std::vector<bool>> diag_matrix;

    std::vector<std::vector<int>> values_matrix;
public:
    NeedlemanWunschAlignment(std::string s, std::string t)
    {
        this->s = s;
        this->t = t;

        up_matrix = std::vector<std::vector<bool>>(s.size() + 1, std::vector<bool>(t.size() + 1, false));
        left_matrix = std::vector<std::vector<bool>>(s.size() + 1, std::vector<bool>(t.size() + 1, false));
        diag_matrix = std::vector<std::vector<bool>>(s.size() + 1, std::vector<bool>(t.size() + 1, false));
        values_matrix = std::vector<std::vector<int>>(s.size() + 1, std::vector<int>(t.size() + 1, 0));

        for (int i = 0; i < s.size() + 1; i++)
        {
            values_matrix[i][0] = i * -2;
            up_matrix[i][0] = true;   
        }
        for (int i = 0; i < t.size() + 1; i++)
        {
            values_matrix[0][i] = i * -2;
            left_matrix[0][i] = true;   
        }

        int rows = s.size();
        int cols = t.size();

        for (int i = 1; i < s.size() + 1; i++) 
        {
            for (int j = 1; j < t.size() + 1; j++)
            {
                int new_val;
                if (s[i - 1] == t[j - 1])
                {
                    new_val = std::max({values_matrix[i - 1][j] - 2, values_matrix[i][j - 1] - 2, values_matrix[i - 1][j - 1] + 1});
                }
                else
                {
                    new_val = std::max({values_matrix[i - 1][j] - 2, values_matrix[i][j - 1] - 2, values_matrix[i - 1][j - 1] - 1});
                }
                if (new_val == values_matrix[i - 1][j] - 2)
                {
                    up_matrix[i][j] = true;
                }
                if (new_val == values_matrix[i][j - 1] - 2)
                {
                    left_matrix[i][j] = true;
                }
                if ((new_val == values_matrix[i - 1][j - 1] + 1 && s[i - 1] == t[j - 1]) || (new_val  == values_matrix[i - 1][j - 1] - 1 && s[i - 1] != t[j - 1]))
                {
                    diag_matrix[i][j] = true;
                }

                values_matrix[i][j] = new_val;
            }
        }

    }

    std::vector<std::vector<int>> get_score_matrix()
    {
        return this->values_matrix;
    }

    int get_score()
    {
        return this->values_matrix[s.size()][t.size()];
    }

    std::vector<std::pair<std::string, std::string>> get_allignments()
    {
        return this->get_allignments(this->s.size(), this->t.size());
    }

    unsigned long long optimal_allignments_num()
    {
        std::vector<std::vector<unsigned long long>> dp(s.size() + 1, std::vector<unsigned long long>(t.size() + 1, 0));
        int rows = s.size();
        int cols = t.size();

        for (int diagonal = 0; diagonal < rows + cols - 1; diagonal++)
        {
        
            int start_row = std::max(0, diagonal - cols + 1);
            int start_col = std::min(diagonal, cols - 1);

            int diag_len = std::min(rows - start_row - 1, start_col) + 1;

            for (int i = 0; i < s.size() + 1; i++)
            {
                dp[i][0] = 1;
            }
            for (int i = 0; i < t.size() + 1; i++)
            {
                dp[0][i] = 1;
            }

            for (int k = 0; k < diag_len; k++) 
            {
                int i = start_row + 1 + k, j = start_col + 1 - k;
                if (up_matrix[i][j])
                {
                    dp[i][j] += dp[i - 1][j];
                }
                if (left_matrix[i][j])
                {
                    dp[i][j] += dp[i][j - 1];
                }
                if (diag_matrix[i][j])
                {
                    dp[i][j] += dp[i - 1][j - 1];
                }
            }
        }

        return dp[rows][cols];
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

    return std::make_pair(s, t); 
}

int main()
{
    srand(time(NULL));
    std::string s = get_random_dna_sequence(5000);
    std::string t = get_random_dna_sequence(5000);

    std::ofstream file("result_seq.txt");
    const auto start{std::chrono::steady_clock::now()};
    NeedlemanWunschAlignment alineacion(s, t);
    const auto end{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_seconds{end - start};
    std::cout << elapsed_seconds.count() << '\n';

    // std::vector<std::vector<int>> score_matrix = alineacion.get_score_matrix();
    // std::vector<std::pair<std::string, std::string>> found_alignments = alineacion.get_allignments();


    file << "Resultados para las cadenas:\n";
    file << "S: " << s << '\n';
    file << "T: " << t << "\n\n";
    file << "Score: " << alineacion.get_score() << '\n';
    // file << "Matriz de Scores\n";

    // file << " \t-\t\t";
    // for (int i = 0; i < t.size(); i++)
    // {
    //     file << t[i] << "\t\t";
    // }

    // file << '\n';
    // for (int i = 0; i < score_matrix.size(); i++)
    // {
    //      std::cout << "que\n";

    //     if (i == 0) file << "-";
    //     else file << s[i - 1];
    //     for (int j = 0; j < score_matrix[i].size(); j++)
    //     {
    //         if (alineacion.diag_matrix[i][j]) file << "\\\t";
    //         else file << "\t";
    //         if (alineacion.up_matrix[i][j]) file << "^\t";
    //         else file << "\t";
    //     }
    //     file << '\n';

    //     for (int j = 0; j < score_matrix[i].size(); j++)
    //     {
    //         if (alineacion.left_matrix[i][j]) file << "<\t";
    //         else file << "\t";
    //         file << score_matrix[i][j] << '\t';
    //     }
    //     file << '\n';
    // }
    // file << '\n';

    // file << "Cantidad de alineamientos producidos: " << alineacion.optimal_allignments_num() << '\n';
    // file << "Alineamientos generados:\n";


    // for (int i = 0; i < found_alignments.size(); i++)
    // {
        // file << found_alignments[i].first << '\n' << found_alignments[i].second << "\n\n";
    // }

    file.close();


    return 0;
}