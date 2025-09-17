// Algotester Library
// Version 1.0.497

#pragma once
#include <iostream>
#include <climits>
#include <stdexcept>
#include <cctype>
#include <regex>
#include <assert.h>
#include <stdarg.h>
#include <array>
#include <tuple>
#include <fstream>
#include <csignal>
#include <cmath>

using namespace std;

#ifdef ALGOTESTER_TEST
#define EXIT(msg) throw runtime_error(msg)
#else
#define EXIT(msg) cout << msg << endl; exit(0)
#endif

class AlgotesterReader;

/// @brief Initialize validator.
/// @return Reader that should be used to read test data.
AlgotesterReader initValidator();

/// @brief Initialize checker.
/// @return Array of readers in order: test input reader, user output reader, correct answer reader.
array<AlgotesterReader, 3> initChecker(int argc, char* argv[]);

/// @brief Initialize interactive checker.
/// @return Tuple of test input reader, user input stream and `AlgotesterReader` from user output stream.
tuple<AlgotesterReader, ofstream, AlgotesterReader> initInteractiveChecker(int argc, char* argv[]);

/// @brief Initialize scorer.
/// @return Array of readers in order: test input reader, user output reader, correct answer reader.
array<AlgotesterReader, 3> initScorer(int argc, char* argv[]);

/// @brief Check that 2 double values are close with given absolute and relative tolerance.
/// @param actual Actual value.
/// @param expected Expected value.
/// @param atol Absolute tolerance.
/// @param rtol Relative tolerance.
/// @return Whether given values are close.
bool isClose(long double actual, long double expected, long double atol=1e-7, long double rtol=1e-7);

/// @brief Case insensitive compare of 2 strings.
/// @param x First string.
/// @param y Second string.
/// @return Whether given strings are equal when case is ignored.
bool equalsIgnoreCase(string x, string y);

/// @brief Check if condition holds and exit with message if not.
void check(bool condition, string message);

class AlgotesterReader{
public:
    AlgotesterReader(istream* _in=&cin, bool _strict=true): in(_in), strict(_strict), isSavedChar(false), lastChar(-1), row(1), col(1) {}
    AlgotesterReader(string filename, bool _strict=true): strict(_strict), isSavedChar(false), lastChar(-1), row(1), col(1) {
        in = new ifstream(filename);
    }

    /// @brief Read single char and check that its value is among possible values.
    /// @param values If non-empty, string that contain all possible char values.
    /// @param skipWhitespace Whether to skip whitespaces before reading next character.
    /// @param name Name of the character.
    /// @return Value that was read.
    /// @note `skipWhitespace=true` could not be used if `strict=true`
    char readChar(string values="", bool skipWhitespace=false, string name="");

    /// @brief Read space character.
    void readSpace() { readChar(" "); }

    /// @brief Read endline character.
    void readEndl();

    /// @brief Read end of file character.
    void readEof();

    /// @brief Check if it is end of file.
    /// @param skipWhitespace If true, whitespaces are skipped, i.e. result of this function is true if only whitespaces are left.
    /// @return If the reader reached end of file.
    /// @note `skipWhitespace=true` could not be used if `strict=true`
    bool isEof(bool skipWhitespace=false);
    
    /// @brief Read long long.
    /// @param name Name of the variable.
    /// @return Value that was read.
    long long readInt(string name="");

    /// @brief Read long long and check if it's in given range.
    /// @param l Minimal value of integer.
    /// @param r Maximal value of integer.
    /// @param name Name of the variable.
    /// @return Value that was read.
    long long readInt(long long l, long long r, string name="");

    /// @brief Read vector of long long with given length.
    /// @param n Length of vector.
    /// @param name Name 
    /// @return Values that were read.
    /// @note In strict mode elements should be space-separated with endl in the end.
    vector<long long> readVectorInt(int n, string name="");

    /// @brief Read vector of long long with given length and check that each element is in given range.
    /// @param n Length of vector.
    /// @param l Minimal value of element.
    /// @param r Maximal value of element.
    /// @param name Name 
    /// @return Values that were read.
    /// @note In strict mode elements should be space-separated with endl in the end.
    vector<long long> readVectorInt(int n, long long l, long long r, string name="");

    /// @brief Read unsigned long long.
    /// @param name Name of the variable.
    /// @return Value that was read.
    unsigned long long readUInt(string name="");

    /// @brief Read unsigned long long and check if it's in given range.
    /// @param l Minimal value.
    /// @param r Maximal value.
    /// @param name Name of the variable.
    /// @return Value that was read.
    unsigned long long readUInt(unsigned long long l, unsigned long long r, string name="");

    /// @brief Read token till the next whitespace and check if it matches given regex.
    /// @param name Name of the variable.
    /// @param regex Regular expression that token must match.
    /// @return Token that was read.
    string readToken(string name="", string regex=".*");

    /// @brief Read token till the next whitespace and check if it matches given regex and its length is in given range.
    /// @param minLength Minimal length of token.
    /// @param maxLength Maximal length of token.
    /// @param name Name of the variable.
    /// @param regex Regular expression that token must match.
    /// @return Token that was read.
    string readToken(int minLength, int maxLength, string name="", string regex=".*");

    /// @brief Read line till the next newline character.
    /// @return Line that was read.
    string readLine();

    /// @brief Read long double.
    /// @param name Name of the variable.
    /// @return Value that was read.
    long double readDouble(string name="");

    /// @brief Read long double and check if it's in given range.
    /// @param l Minimal value.
    /// @param r Maximal value.
    /// @param name Name of the variable.
    /// @return Value that was read.
    long double readDouble(long double l, long double r, string name="");

    /// @brief Read long double and check its characteristics.
    /// @param minDecimals Minimal number of decimal digits after the point.
    /// @param maxDecimals Maximal number of decimal digits after the point.
    /// @param allowExpFormat Whether exponential format is allowed.
    /// @param name Name of the variable.
    /// @return Value that was read.
    long double readDouble(int minDecimals, int maxDecimals, bool allowExpFormat, string name);

    /// @brief Read long double and check its characteristics and range.
    /// @param l Minimal value.
    /// @param r Maximal value.
    /// @param minDecimals Minimal number of decimal digits after the point.
    /// @param maxDecimals Maximal number of decimal digits after the point.
    /// @param allowExpFormat Whether exponential format is allowed.
    /// @param name Name of the variable.
    /// @return Value that was read.
    long double readDouble(long double l, long double r, int minDecimals, int maxDecimals, bool allowExpFormat, string name);

    /// @brief Set the callback function that is executed before the error message is printed to stdout.
    /// @param error_callback Callback function.
    void setErrorCallback(std::function<void()> error_callback);

private:
    istream* in;
    bool strict;
    bool isSavedChar;
    char lastChar;
    int row, col;
    int row_prev, col_prev;
    std::function<void()> error_callback = [](){};

    char getChar() {
        char res;
        if (isSavedChar) {
            res = lastChar;
            isSavedChar = false;
        } else {
            if (!in->get(res)) {
                return 0;
            }
        }
        row_prev = row;
        col_prev = col;
        if (res == '\n') {
            row ++;
            col = 1;
        } else {
            col ++;
        }
        return res;
    }

    void error(string message) {
        error_callback();
        EXIT(message);
    }

    string error_message(const char* s, ...) {
        va_list arg;
        char res[128];
        va_start(arg, s);
        vsnprintf (res, 128, s, arg);
        va_end(arg);
        return this->info() + string(res);
    }

    string info() { return "row " + std::to_string(row) + ": col " + std::to_string(col) + ": "; }

    string readTokenImpl(int minLength, int maxLength, string name="") {
        string res = "";
        if (!strict)
            readWhitespaces();
        while (true) {
            char ch = getChar();
            if (!ch || isspace(ch)) {
                if (ch)
                    saveChar(ch);
                if (res.empty()) {
                    error(error_message("Token %s is empty (maybe you forget to read whitespace character)", name.c_str()));
                }
                if ((int)res.size() < minLength) {
                    error(error_message("Token %s is too short, expected length in range [%d, %d]", name.c_str(), minLength, maxLength));
                }
                return res;
            }
            if ((int)res.size() >= maxLength) {
                error(error_message("Token %s is too long, expected length in range [%d, %d]", name.c_str(), minLength, maxLength));
            }
            res += ch;
        }
    }

    void match_string(string s, string r) {
        if (r == ".*") {
            return;
        }
        const std::regex regex(r);
        if (!std::regex_match(s, regex)) {
            error(error_message("Token %s doesn't match regex %s", s.c_str(), r.c_str()));
        }
    }

    void saveChar(char c) {
        row = row_prev;
        col = col_prev;
        assert(!isSavedChar);
        lastChar = c;
        isSavedChar = true;
    }

    unsigned long long parseUInt(string s, unsigned long long limit, bool neg, string name = "") {
        string fs = (neg ? "-" : "") + s;
        if (!strict) {
            int l = 0;
            while (s.size() - l > 1 && s[l] == '0') ++ l;
            s = s.substr(l);
        }
        if (s.size() > 20) {
            error(error_message("Integer %s out of long long range: %s", name.c_str(), fs.c_str()));
        }
        unsigned long long res = 0;
        for(unsigned i = 0; i < s.size(); ++i) {
            int digit = s[i] - '0';
            if (digit < 0 || digit > 9) {
                error(error_message("Cannot parse %s: %s", name.c_str(), fs.c_str()));
            }
            if (strict && i == 0 && digit == 0 && (s.size() > 1 || neg)) {
                if (neg && s.size() == 1) {
                    error(error_message("Negative zero is not allowed in integer %s", name.c_str()));
                }
                error(error_message("Leading zeros in integer %s: %s", name.c_str(), fs.c_str()));
            }
            if ((limit - digit) / 10 < res) {
                error(error_message("Integer %s out of long long range: %s", name.c_str(), fs.c_str()));
            }
            res = res * 10 + digit;
        }
        return res;
    }

    void readWhitespaces() {
        while(true) {
            char ch = getChar();
            if (ch == 0)
                break;
            if (isspace(ch))
                continue;
            saveChar(ch);
            break;
        }
    }

    string charInfo(char ch) {
        string res = "";
        switch (ch)
        {
        case '\n':
            res += "\\n";
            break;
        case '\r':
            res += "\\r";
            break;
        case ' ':
            res += "space";
            break;
        default:
            res += ch;
            break;
        }
        res += "(" + to_string(int(ch)) + ")";
        return res;
    }
};

void check(bool condition, string message, std::function<void()> error_callback) {
    if (!condition) {
        error_callback();
        EXIT(message);
    }
}

void check(bool condition, string message) {
    check(condition, message, [](){});
}

void AlgotesterReader::setErrorCallback(std::function<void()> error_callback) {
    this->error_callback = error_callback;
}

char AlgotesterReader::readChar(string values, bool skipWhitespace, string name) {
    assert (!strict || !skipWhitespace);
    if (skipWhitespace)
        readWhitespaces();
    char ch = getChar();
    if (ch == 0) {
        error(error_message("Character %s could not be read", name.c_str()));
    }
    if (values.size() != 0) {
        bool ok = false;
        for (char t: values)
            if (ch == t) {
                ok = true;
                break;
            }
        if (!ok) {
            string expected_info = "";
            for(char c: values) {
                expected_info += " " + charInfo(c);
            }
            error(error_message("Wrong character %s, expected values =%s, actual=%s", name.c_str(), expected_info.c_str(), charInfo(ch).c_str()));
        }
    }
    return ch;
}

void AlgotesterReader::readEndl() {
    if (strict) {
        readChar("\n");
    } else {
        char c0;
        while (true) {
            c0 = readChar(" \n\r");
            if (c0 != ' ')
                break;
        } 
        char c1 = getChar();
        if ((c1 != '\n' && c1 != '\r' && c1) || c0 == c1) {
            saveChar(c1);
        }
    }
}

void AlgotesterReader::readEof() {
    if (!strict)
        readWhitespaces();
    char ch = getChar();
    if (ch == 0)
        return;
    error(error_message("Not end of file, read character %s", charInfo(ch).c_str()));
};

bool AlgotesterReader::isEof(bool skipWhitespace) {
    assert(!strict || !skipWhitespace); // Could not use `skipWhitespace=true` in strict mode.
    if (skipWhitespace)
        readWhitespaces();
    char ch = getChar();
    if (ch == 0)
        return true;
    saveChar(ch);
    return false;
}

long long AlgotesterReader::readInt(string name) {
    string s = readToken(name);
    long long res;
    if (s[0] == '-') {
        res = -(long long)parseUInt(s.substr(1), (unsigned long long)LLONG_MAX + 1, true, name);
    } else {
        if (!strict && s[0] == '+') {
            s = s.substr(1);
        }
        res = parseUInt(s, LLONG_MAX, false, name);
    }
    return res;
}

long long AlgotesterReader::readInt(long long l, long long r, string name){
    long long res = readInt(name);
    if (res < l || res > r) {
        error(error_message("Integer %s out of range, expected range=[%lld, %lld], actual=%lld\n", name.c_str(), l, r, res));
    }
    return res;
}

vector<long long> AlgotesterReader::readVectorInt(int n, string name) {
    vector<long long> res;
    for(int i = 0; i < n; ++i) {
        res.push_back(readInt(name + "[" + to_string(i) + "]"));
        if (strict) {
            if (i + 1 < n)
                readSpace();
            else
                readEndl();
        }
    }
    return res;
}

vector<long long> AlgotesterReader::readVectorInt(int n, long long l, long long r, string name) {
    vector<long long> res;
    for(int i = 0; i < n; ++i) {
        res.push_back(readInt(l, r, name + "[" + to_string(i) + "]"));
        if (strict) {
            if (i + 1 < n)
                readSpace();
            else
                readEndl();
        }
    }
    return res;
}

unsigned long long AlgotesterReader::readUInt(string name) {
    string s = readToken();
    unsigned long long res = parseUInt(s, ULLONG_MAX, false, name);
    return res;
}

unsigned long long AlgotesterReader::readUInt(unsigned long long l, unsigned long long r, string name){
    unsigned long long res = this->readUInt(name);
    if (res < l || res > r) {
        error(error_message("Integer %s out of range, expected range=[%llu, %llu], actual=%llu\n", name.c_str(), l, r, res));
    }
    return res;
}

string AlgotesterReader::readToken(string name, string regex) {
    return readToken(0, INT_MAX, name, regex);
}

string AlgotesterReader::readToken(int minLength, int maxLength, string name, string regex) {
    string res = readTokenImpl(minLength, maxLength, name);
    match_string(res, regex);
    return res;
}

string AlgotesterReader::readLine() {
    string res = "";
    while (true) {
        char ch = getChar();
        if (!ch || ch == '\n' || ch == '\r') {
            saveChar(ch);
            if (ch != 0)
                readEndl();
            return res;
        }
        res += ch;
    }
}

long double AlgotesterReader::readDouble(string name) {
    return readDouble(0, 47, !strict, name);
}

long double AlgotesterReader::readDouble(long double l, long double r, string name) {
    long double res = readDouble(name);
    if (res < l || res > r) {
        error(error_message("Double %s out of range, expected range=[%llf, %llf], actual=%llf\n", name.c_str(), l, r, res));
    }
    return res;
}

long double AlgotesterReader::readDouble(int minDecimals, int maxDecimals, bool allowExpFormat, string name) {
    string s = readToken(name);
    long double res;
    try{
        size_t pos;
        res = stold(s, &pos);
        if (pos != s.size()) {
            throw invalid_argument("");
        }
    }
    catch (invalid_argument&) {
        error(error_message("Cannot parse double %s: %s", name.c_str(), s.c_str()));
    }
    catch(out_of_range&) {
        error(error_message("Double %s is out of long double range: %s", name.c_str(), s.c_str()));
    }
    if (!isfinite(res)){
        error(error_message("Double %s is not finite: %s", name.c_str(), s.c_str()));
    }
    int expPos = max((int)s.find('e'), (int)s.find('E'));
    int dotPos = (int)s.find('.');
    if (!allowExpFormat && expPos != -1) {
        error(error_message("Double %s is in exponential format: %s", name.c_str(), s.c_str()));
    }
    int decimals = (dotPos == -1 ? 0 : (expPos == -1 ? (int)s.size() - dotPos - 1 : expPos - dotPos - 1));
    if (decimals < minDecimals || decimals > maxDecimals) {
        error(error_message("Double %s has incorrect number of decimals, expected range=[%d, %d], actual=%d: %s", name.c_str(), minDecimals, maxDecimals, decimals, s.c_str()));
    }
    return res;
}

long double AlgotesterReader::readDouble(long double l, long double r, int minDecimals, int maxDecimals, bool allowExpFormat, string name) {
    long double res = readDouble(minDecimals, maxDecimals, allowExpFormat, name);
    if (res < l || res > r) {
        error(error_message("Double %s out of range, expected range=[%llf, %llf], actual=%llf\n", name.c_str(), l, r, res));
    }
    return res;
}

AlgotesterReader initValidator() {
    return AlgotesterReader();
}

array<AlgotesterReader, 3> initChecker(int argc, char* argv[]) {
    static_cast<void>(argc);
    return {AlgotesterReader(argv[1], false), AlgotesterReader(argv[2], false), AlgotesterReader(argc >= 4 ? argv[3] : "/dev/null", false)};
}

#if defined(__unix__)
#define IGNORE_BROKEN_PIPE signal(SIGPIPE, SIG_IGN)
#else
#define IGNORE_BROKEN_PIPE
#endif

tuple<AlgotesterReader, ofstream, AlgotesterReader> initInteractiveChecker(int argc, char* argv[]) {
    IGNORE_BROKEN_PIPE;
    static_cast<void>(argc);
    return {AlgotesterReader(&cin, false), ofstream(argv[1]), AlgotesterReader(argv[2], false)};
}

array<AlgotesterReader, 3> initScorer(int argc, char* argv[]) {
    return initChecker(argc, argv);
}

bool isClose(long double actual, long double expected, long double atol, long double rtol) {
    return abs(actual - expected) <= max(rtol * abs(expected), atol);
}

bool equalsIgnoreCase(string x, string y) {
    if (x.size() != y.size())
        return false;
    for (unsigned i = 0; i < x.size(); ++i) {
        if (tolower(x[i]) != tolower(y[i]))
            return false;
    }
    return true;
}

