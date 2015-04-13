#include "bio.hpp"

using namespace std;
namespace Bio {

// this returns a vector filled with the first n_max_tokens
std::vector<std::string> tokenize_partial(const std::string &str,
                                          size_t n_max_tokens,
                                          const std::string &search) {
  std::vector<std::string> tokens;
  tokens.reserve(n_max_tokens);
  std::string::size_type p_last = str.find_first_not_of(search, 0);
  std::string::size_type p_curr = str.find_first_of(search, p_last);

  if (n_max_tokens > 1)
    while ((std::string::npos != p_curr || std::string::npos != p_last) &&
           tokens.size() + 1 != n_max_tokens) {
      tokens.push_back(str.substr(p_last, p_curr - p_last));
      p_last = str.find_first_not_of(search, p_curr);
      p_curr = str.find_first_of(search, p_last);
    }

  tokens.push_back(str.substr(p_last, p_curr - p_last));
  return tokens;
}
}
