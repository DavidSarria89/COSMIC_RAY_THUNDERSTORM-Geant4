

#include <G4ThreeVector.hh>
#include "myUtils.hh"

/////////////////////////////////////////////////////////

namespace myUtils {
    ///
    // output in microseconds
    double get_wall_time() {
        std::chrono::high_resolution_clock m_clock;
        double time = std::chrono::duration_cast<std::chrono::microseconds>(
                m_clock.now().time_since_epoch())
                .count();
        return time;
    }

/////////////////////////////////////////////////////////

    ///
    long reversDigits(long num) {
        long rev_num = 0;
        while (num > 0) {
            rev_num = rev_num * 10 + num % 10;
            num = num / 10;
        }
        return rev_num;
    }

    bool addOvf(long *result, long a, long b) {
        // adding two long integers but checking if the sum does not overflow
        if (a > std::numeric_limits<long>::max() - b)
            return false;
        else {
            *result = a + b;
            return true;
        }
    }

    std::string convertToString(char a[37], int size) {
        int i;
        std::string s = "";
        for (i = 0; i < size; i++) {
            s = s + a[i];
        }
        return s;
    }

    std::string change_letters_to_numbers(std::string my_str) {
        std::replace(my_str.begin(), my_str.end(), 'a', '1');
        std::replace(my_str.begin(), my_str.end(), 'b', '2');
        std::replace(my_str.begin(), my_str.end(), 'c', '3');
        std::replace(my_str.begin(), my_str.end(), 'd', '4');
        std::replace(my_str.begin(), my_str.end(), 'e', '5');
        std::replace(my_str.begin(), my_str.end(), 'f', '6');
        std::replace(my_str.begin(), my_str.end(), 'g', '7');
        std::replace(my_str.begin(), my_str.end(), 'h', '8');
        std::replace(my_str.begin(), my_str.end(), 'i', '9');
        std::replace(my_str.begin(), my_str.end(), 'j', '1');
        std::replace(my_str.begin(), my_str.end(), 'j', '2');
        std::replace(my_str.begin(), my_str.end(), 'l', '3');
        std::replace(my_str.begin(), my_str.end(), 'm', '4');
        std::replace(my_str.begin(), my_str.end(), 'n', '5');
        std::replace(my_str.begin(), my_str.end(), 'o', '6');
        std::replace(my_str.begin(), my_str.end(), 'p', '7');
        std::replace(my_str.begin(), my_str.end(), 'k', '8');
        std::replace(my_str.begin(), my_str.end(), 'r', '9');
        std::replace(my_str.begin(), my_str.end(), 's', '1');
        std::replace(my_str.begin(), my_str.end(), 't', '2');
        std::replace(my_str.begin(), my_str.end(), 'u', '3');
        std::replace(my_str.begin(), my_str.end(), 'v', '4');
        std::replace(my_str.begin(), my_str.end(), 'w', '5');
        std::replace(my_str.begin(), my_str.end(), 'x', '6');
        std::replace(my_str.begin(), my_str.end(), 'y', '7');
        std::replace(my_str.begin(), my_str.end(), 'z', '8');
        my_str.erase(remove(my_str.begin(), my_str.end(), '-'), my_str.end());
        return my_str;
    }

    long generate_a_unique_ID() {

        uuid_t uu;
        uuid_generate(uu);
        char uuid[37];
        uuid_unparse(uu, uuid);
//        std::cout << uuid << std::endl;
        int mys = sizeof(uuid) / sizeof(char);
        std::string my_str = std::string(convertToString(uuid, mys));

        my_str = change_letters_to_numbers(my_str);

        my_str = my_str.substr(1, my_str.size() - 24);
        long output = std::stol(my_str);
//        int dummy = 1 + 1;
        std::cout << output << std::endl;
        return output;
    }


    int generate_a_unique_ID_int32() {

        uuid_t uu;
        uuid_generate(uu);
        char uuid[37];
        uuid_unparse(uu, uuid);
//        std::cout << uuid << std::endl;
        int mys = sizeof(uuid) / sizeof(char);
        std::string my_str = std::string(convertToString(uuid, mys));

        my_str = change_letters_to_numbers(my_str);

        my_str = my_str.substr(1, my_str.size() - 25);
        long output = std::stol(my_str);
//        int dummy = 1 + 1;
        std::cout << output << std::endl;

        if (output>INT_MAX) std::abort();

        return output;
    }

    ///

    std::vector<double> logspace(const double a, const double b, const int n) {
        double c;
        int i;

        std::vector<double> output;
        output.clear();
        output.reserve(n);

        /* step size */
        c = (b - a) / (n - 1);

        /* fill vector */
        for (i = 0; i < n - 1; ++i)
            output.push_back(pow(10., a + i * c));

        /* fix last entry to 10^b */
        output.push_back(pow(10., b));

        return output;
    }

    double check_available_RAM() { // output in GB
        struct sysinfo memInfo;

        sysinfo (&memInfo);
        long long totalVirtualMem = memInfo.totalram;
//Add other values in next statement to avoid int overflow on right hand side...
        totalVirtualMem += memInfo.totalswap;
        totalVirtualMem *= memInfo.mem_unit;

        long long totalPhysMem = memInfo.totalram;
//Multiply in next statement to avoid int overflow on right hand side...
        totalPhysMem *= memInfo.mem_unit;

        return double(totalPhysMem) / 1000 / 1000 / 1000;
    }

    double check_USED_RAM() {

        FILE *file = fopen("/proc/self/status", "r");
        int result = -1;
        char line[128];

        while (fgets(line, 128, file) != NULL) {
            if (strncmp(line, "VmRSS:", 6) == 0) {
                result = parseLine(line);
                break;
            }
        }
        fclose(file);
        const double arbitrary_factor = 0.82;

        return double(result) / 1024 / 1024 * arbitrary_factor;
    }

    int parseLine(char *line) {
        // This assumes that a digit will be found and the line ends in " Kb".
        int i = strlen(line);
        const char *p = line;
        while (*p < '0' || *p > '9') p++;
        line[i - 3] = '\0';
        i = atoi(p);
        return i;
    }


}

