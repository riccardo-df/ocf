#include <math.h>
#include <iostream>
#include <sstream>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "utility.h"
#include "globals.h"
#include "Data.h"

namespace ocf {

void equalSplit(std::vector<unsigned int>& result, unsigned int start, unsigned int end, unsigned int num_parts) {

  result.reserve(num_parts + 1);

  // Return range if only 1 part.
  if (num_parts == 1) {
    result.push_back(start);
    result.push_back(end + 1);
    
    return;
  }

  // Return vector from start to end+1 if more parts than elements.
  if (num_parts > end - start + 1) {
    for (unsigned int i = start; i <= end + 1; ++i) {
      result.push_back(i);
    }
    
    return;
  }

  unsigned int length = (end - start + 1);
  unsigned int part_length_short = length / num_parts;
  unsigned int part_length_long = (unsigned int) ceil(length / ((double) num_parts));
  unsigned int cut_pos = length % num_parts;

  // Adding long ranges.
  for (unsigned int i = start; i < start + cut_pos * part_length_long; i = i + part_length_long) {
    result.push_back(i);
  }

  // Adding short ranges.
  for (unsigned int i = start + cut_pos * part_length_long; i <= end + 1; i = i + part_length_short) {
    result.push_back(i);
  }
}

void loadDoubleVectorFromFile(std::vector<double>& result, std::string filename) { // #nocov start
  // Opening input file.
  std::ifstream input_file;
  input_file.open(filename);
  
  if (!input_file.good()) {
    throw std::runtime_error("Could not open file: " + filename);
  }

  // Reading the first line, ignoring the rest.
  std::string line;
  getline(input_file, line);
  std::stringstream line_stream(line);
  double token;
  
  while (line_stream >> token) {
    result.push_back(token);
  }
} // #nocov end

void drawWithoutReplacement(std::vector<size_t>& result, std::mt19937_64& random_number_generator, size_t max,
    size_t num_samples) {
  if (num_samples < max / 10) {
    drawWithoutReplacementSimple(result, random_number_generator, max, num_samples);
  } else {
    drawWithoutReplacementFisherYates(result, random_number_generator, max, num_samples);
  }
}

void drawWithoutReplacementSkip(std::vector<size_t>& result, std::mt19937_64& random_number_generator, size_t max,
    const std::vector<size_t>& skip, size_t num_samples) {
  if (num_samples < max / 10) {
    drawWithoutReplacementSimple(result, random_number_generator, max, skip, num_samples);
  } else {
    drawWithoutReplacementFisherYates(result, random_number_generator, max, skip, num_samples);
  }
}

void drawWithoutReplacementSimple(std::vector<size_t>& result, std::mt19937_64& random_number_generator, size_t max,
    size_t num_samples) {

  result.reserve(num_samples);

  // Setting all to not selected.
  std::vector<bool> temp;
  temp.resize(max, false);

  std::uniform_int_distribution<size_t> unif_dist(0, max - 1);
  
  for (size_t i = 0; i < num_samples; ++i) {
    size_t draw;
  
    do {
      draw = unif_dist(random_number_generator);
    } while (temp[draw]);
    temp[draw] = true;
    result.push_back(draw);
  }
}

void drawWithoutReplacementSimple(std::vector<size_t>& result, std::mt19937_64& random_number_generator, size_t max,
    const std::vector<size_t>& skip, size_t num_samples) {

  result.reserve(num_samples);

  // Setting all to not selected.
  std::vector<bool> temp;
  temp.resize(max, false);

  std::uniform_int_distribution<size_t> unif_dist(0, max - 1 - skip.size());
  
  for (size_t i = 0; i < num_samples; ++i) {
    size_t draw;
    
    do {
      draw = unif_dist(random_number_generator);
      
      for (auto& skip_value : skip) {
        if (draw >= skip_value) {
          ++draw;
        }
      }
    } while (temp[draw]);
    temp[draw] = true;
    result.push_back(draw);
  }
}

void drawWithoutReplacementFisherYates(std::vector<size_t>& result, std::mt19937_64& random_number_generator,
    size_t max, size_t num_samples) {

  // Creating indices.
  result.resize(max);
  std::iota(result.begin(), result.end(), 0);

  // Drawing without replacement using Fisher Yates algorithm.
  std::uniform_real_distribution<double> distribution(0.0, 1.0);
  for (size_t i = 0; i < num_samples; ++i) {
    size_t j = i + distribution(random_number_generator) * (max - i);
    std::swap(result[i], result[j]);
  }

  result.resize(num_samples);
}

void drawWithoutReplacementFisherYates(std::vector<size_t>& result, std::mt19937_64& random_number_generator,
    size_t max, const std::vector<size_t>& skip, size_t num_samples) {

  // Creating indices.
  result.resize(max);
  std::iota(result.begin(), result.end(), 0);

  // Skipping indices.
  for (size_t i = 0; i < skip.size(); ++i) {
    result.erase(result.begin() + skip[skip.size() - 1 - i]);
  }

  // Drawing without replacement using Fisher Yates algorithm.
  std::uniform_real_distribution<double> distribution(0.0, 1.0);
  for (size_t i = 0; i < num_samples; ++i) {
    size_t j = i + distribution(random_number_generator) * (max - skip.size() - i);
    std::swap(result[i], result[j]);
  }

  result.resize(num_samples);
}

void drawWithoutReplacementWeighted(std::vector<size_t>& result, std::mt19937_64& random_number_generator,
    size_t max_index, size_t num_samples, const std::vector<double>& weights) {

  result.reserve(num_samples);

  // Setting all to not selected.
  std::vector<bool> temp;
  temp.resize(max_index + 1, false);

  std::discrete_distribution<> weighted_dist(weights.begin(), weights.end());
  
  for (size_t i = 0; i < num_samples; ++i) {
    size_t draw;
    
    do {
      draw = weighted_dist(random_number_generator);
    } while (temp[draw]);
    temp[draw] = true;
    result.push_back(draw);
  }
}

std::string uintToString(unsigned int number) {
#if WIN_R_BUILD == 1
  std::stringstream temp;
  temp << number;
  return temp.str();
#else
  return std::to_string(number);
#endif
}

std::string beautifyTime(unsigned int seconds) { // #nocov start
  std::string result;

  // Adding seconds, minutes, hours, days if larger than zero.
  unsigned int out_seconds = (unsigned int) seconds % 60;
  result = uintToString(out_seconds) + " seconds";
  unsigned int out_minutes = (seconds / 60) % 60;
  
  if (seconds / 60 == 0) {
    return result;
  } else if (out_minutes == 1) {
    result = "1 minute, " + result;
  } else {
    result = uintToString(out_minutes) + " minutes, " + result;
  }
  
  unsigned int out_hours = (seconds / 3600) % 24;
  
  if (seconds / 3600 == 0) {
    return result;
  } else if (out_hours == 1) {
    result = "1 hour, " + result;
  } else {
    result = uintToString(out_hours) + " hours, " + result;
  }
  
  unsigned int out_days = (seconds / 86400);
  
  if (out_days == 0) {
    return result;
  } else if (out_days == 1) {
    result = "1 day, " + result;
  } else {
    result = uintToString(out_days) + " days, " + result;
  }
  return result;
} // #nocov end

// #nocov start
size_t roundToNextMultiple(size_t value, unsigned int multiple) {

  if (multiple == 0) {
    return value;
  }

  size_t remainder = value % multiple;
  if (remainder == 0) {
    return value;
  }

  return value + multiple - remainder;
}
// #nocov end

void splitString(std::vector<std::string>& result, const std::string& input, char split_char) { // #nocov start

  std::istringstream ss(input);
  std::string token;

  while (std::getline(ss, token, split_char)) {
    result.push_back(token);
  }
} // #nocov end

void splitString(std::vector<double>& result, const std::string& input, char split_char) { // #nocov start

  std::istringstream ss(input);
  std::string token;

  while (std::getline(ss, token, split_char)) {
    result.push_back(std::stod(token));
  }
} // #nocov end

void shuffleAndSplit(std::vector<size_t>& first_part, std::vector<size_t>& second_part, size_t n_all, size_t n_first,
    std::mt19937_64 random_number_generator) {

  // Reserving space.
  first_part.resize(n_all);

  // Filling with 0, ..., n_all-1 and shuffing.
  std::iota(first_part.begin(), first_part.end(), 0);
  std::shuffle(first_part.begin(), first_part.end(), random_number_generator);

  // Copying to second part.
  second_part.resize(n_all - n_first);
  std::copy(first_part.begin() + n_first, first_part.end(), second_part.begin());

  // Resizing first part.
  first_part.resize(n_first);
}

void shuffleAndSplitAppend(std::vector<size_t>& first_part, std::vector<size_t>& second_part, size_t n_all,
    size_t n_first, const std::vector<size_t>& mapping, std::mt19937_64 random_number_generator) {
  // Old end is start position for new data.
  size_t first_old_size = first_part.size();
  size_t second_old_size = second_part.size();

  // Reserving space.
  first_part.resize(first_old_size + n_all);
  std::vector<size_t>::iterator first_start_pos = first_part.begin() + first_old_size;

  // Filling with 0, ..., n_all-1 and shuffling.
  std::iota(first_start_pos, first_part.end(), 0);
  std::shuffle(first_start_pos, first_part.end(), random_number_generator);

  // Mapping.
  for (std::vector<size_t>::iterator j = first_start_pos; j != first_part.end(); ++j) {
    *j = mapping[*j];
  }

  // Copying to second part.
  second_part.resize(second_part.size() + n_all - n_first);
  std::vector<size_t>::iterator second_start_pos = second_part.begin() + second_old_size;
  std::copy(first_start_pos + n_first, first_part.end(), second_start_pos);

  // Resizing first part.
  first_part.resize(first_old_size + n_first);
}

std::string checkUnorderedVariables(const Data& data, const std::vector<std::string>& unordered_variable_names) { // #nocov start
  size_t num_rows = data.getNumRows();
  std::vector<size_t> sampleIDs(num_rows);
  std::iota(sampleIDs.begin(), sampleIDs.end(), 0);

  // Checking for all unordered variables.
  for (auto& variable_name : unordered_variable_names) {
    size_t varID = data.getVariableID(variable_name);
    std::vector<double> all_values;
    data.getAllValues(all_values, sampleIDs, varID, 0, sampleIDs.size());

    // Checking level count.
    size_t max_level_count = 8 * sizeof(size_t) - 1;
    if (all_values.size() > max_level_count) {
      return "Too many levels in unordered categorical variable " + variable_name + ". Only "
          + uintToString(max_level_count) + " levels allowed on this system.";
    }

    // Checking positive integers.
    if (!checkPositiveIntegers(all_values)) {
      return "Not all values in unordered categorical variable " + variable_name + " are positive integers.";
    }
  }
  
  return "";
} // #nocov end

bool checkPositiveIntegers(const std::vector<double>& all_values) { // #nocov start
  for (auto& value : all_values) {
    if (value < 1 || !(floor(value) == value)) {
      return false;
    }
  }
  return true;
} // #nocov end

std::vector<size_t> numSamplesLeftOfCutpoint(std::vector<double>& x, const std::vector<size_t>& indices) {
  std::vector<size_t> num_samples_left;
  num_samples_left.reserve(x.size());

  for (size_t i = 0; i < x.size(); ++i) {
    if (i == 0) {
      num_samples_left.push_back(1);
    } else if (x[indices[i]] == x[indices[i - 1]]) {
      ++num_samples_left[num_samples_left.size() - 1];
    } else {
      num_samples_left.push_back(num_samples_left[num_samples_left.size() - 1] + 1);
    }
  }

  return num_samples_left;
}

// #nocov start
std::stringstream& readFromStream(std::stringstream& in, double& token) {
  if (!(in >> token) && (std::fpclassify(token) == FP_SUBNORMAL)) {
    in.clear();
  }
  return in;
}
// #nocov end

} // namespace ocf
