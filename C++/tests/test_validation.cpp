// test_validation.cpp — cross-validate C++ vs Python 96th percentile results
//
// Reads two CSV files and compares the 96th percentile values per category.
// Expected to be run after both Python and C++ simulations have produced output.
//
// Usage:
//   test_validation --mode ecai   (compares ecai_results.csv vs ecai_results_py.csv)
//   test_validation --mode prob   (compares ecai_ashrae241_96th_percentile.csv vs prob_py.csv)

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <sstream>

struct CSVRow {
    std::string category;
    double value;
};

static std::vector<CSVRow> load_csv(const std::string& path, int value_col) {
    std::ifstream in(path);
    if (!in) {
        fprintf(stderr, "Cannot open %s\n", path.c_str());
        exit(1);
    }
    std::vector<CSVRow> rows;
    std::string line;
    std::getline(in, line); // skip header
    while (std::getline(in, line)) {
        std::stringstream ss(line);
        std::string field;
        std::vector<std::string> fields;
        while (std::getline(ss, field, ',')) fields.push_back(field);
        if (fields.size() > static_cast<size_t>(value_col)) {
            try {
                double val = std::stod(fields[value_col]);
                rows.push_back({fields[0], val});
            } catch (...) {
                // Skip non-numeric rows
            }
        }
    }
    return rows;
}

int main(int argc, char* argv[]) {
    std::string mode = "ecai";
    double tolerance = 0.05;  // 5% relative tolerance for most categories
    double high_var_tolerance = 0.15;  // 15% for high-variance categories

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--mode" && i + 1 < argc) mode = argv[++i];
        else if (arg == "--tol" && i + 1 < argc) tolerance = std::atof(argv[++i]);
    }

    std::string cpp_file, py_file;
    int value_col;

    if (mode == "ecai") {
        cpp_file = "ecai_results.csv";
        py_file = "ecai_results_py.csv";
        value_col = 1;  // Percentile_96_Lps
    } else {
        cpp_file = "ecai_ashrae241_96th_percentile.csv";
        py_file = "prob_py.csv";
        value_col = 2;  // P_96th_percentile
    }

    auto cpp_rows = load_csv(cpp_file, value_col);
    auto py_rows = load_csv(py_file, value_col);

    // Build maps for lookup
    std::map<std::string, double> cpp_map, py_map;
    for (const auto& r : cpp_rows) cpp_map[r.category] = r.value;
    for (const auto& r : py_rows) py_map[r.category] = r.value;

    // High-variance categories
    std::vector<std::string> high_var = {"Dwelling", "Exam", "Patient", "Resident", "Sorting"};

    int failures = 0;
    int passes = 0;

    printf("%-15s  %12s  %12s  %10s  %s\n", "Category", "C++", "Python", "RelDiff", "Status");
    printf("%s\n", std::string(70, '-').c_str());

    for (const auto& [cat, py_val] : py_map) {
        auto it = cpp_map.find(cat);
        if (it == cpp_map.end()) {
            printf("%-15s  (missing from C++)\n", cat.c_str());
            failures++;
            continue;
        }
        double cpp_val = it->second;

        if (py_val == 0 && cpp_val == 0) {
            printf("%-15s  %12.4f  %12.4f  %9.1f%%  PASS\n", cat.c_str(), cpp_val, py_val, 0.0);
            passes++;
            continue;
        }

        double rel_diff = std::fabs(cpp_val - py_val) / std::fabs(py_val);

        bool is_high_var = false;
        for (const auto& hv : high_var)
            if (hv == cat) { is_high_var = true; break; }

        double tol = is_high_var ? high_var_tolerance : tolerance;
        bool ok = rel_diff < tol;

        printf("%-15s  %12.4f  %12.4f  %9.1f%%  %s%s\n",
               cat.c_str(), cpp_val, py_val, rel_diff * 100,
               ok ? "PASS" : "FAIL",
               is_high_var ? " (high-var)" : "");

        if (ok) passes++;
        else failures++;
    }

    printf("\n%d/%d categories passed (tolerance: %.0f%% normal, %.0f%% high-variance)\n",
           passes, passes + failures, tolerance * 100, high_var_tolerance * 100);

    return failures == 0 ? 0 : 1;
}