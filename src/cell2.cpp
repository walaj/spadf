#include "Cell2.h"
#include <fstream>
#include <vector>
#include <iostream>
#include <cassert>
#include <cstdio>
#include <cstring>
#include <sstream>
#include <algorithm>



#define MAX_LINE_SIZE 2048

std::vector<double> CellTable::XY() const {

  std::vector<double> out;
  out.reserve(2*num_cells());
  assert(!x.empty());
  assert(!y.empty());
  assert(table.at(x).size());
  assert(table.at(y).size());  
  out.insert(out.end(), table.at(x).begin(), table.at(x).end());
  out.insert(out.end(), table.at(y).begin(), table.at(y).end());  
  
  return out;
  
}

void CellTable::FillMarkerNames(char** array) const {
  int i = 0;
  for (const auto& m : markers) {
    array[i] = new char[m.size() + 1];
    strcpy(array[i], m.c_str());
    i++;
  }
}

void CellTable::FillMetaNames(char** array) const {
  int i = 0;
  for (const auto& m : meta) {
    array[i] = new char[m.size() + 1];
    strcpy(array[i], m.c_str());
    i++;    
  }
}

CellTable::CellTable(const char* file, const char* markers_file) {

  // read the markers first
  std::ifstream       mfile(markers_file);
  std::string         line;

  // read in the markers
  std::string mm = std::string(markers_file);
  JsonReader json_reader(mm);
  json_reader.ReadData();
  
  //debug
  std::cerr << " JSON X " << json_reader.GetX() << " JSON Y " << json_reader.GetY() << std::endl;
  
  x = json_reader.GetX();
  y = json_reader.GetY();
  assert(!x.empty());
  assert(!y.empty());  

  markers = json_reader.GetMarkers();
  meta = json_reader.GetMetaCells();  
  
  //std::getline(mfile, line); 
  //std::istringstream stream(line);
  //std::string item;
  //while (std::getline(mfile, item)) {
  //  arkers.insert(item);
  //}
  
  // read in the header
  std::ifstream       rfile(file);  
  std::getline(rfile, line); 
  std::istringstream header_stream(line);
  std::string header_item;

  // store the markers, meta etc in an ordered way
  // this is needed to link the vec's to the right
  // table element, since unordered_map doesn't store
  // by order of insertion
  std::vector<std::string> elements; 
  
  while (std::getline(header_stream, header_item, ',')) {
    elements.push_back(header_item);
    table[header_item] = std::vector<float>();
    //    std::cerr << " H " << header_item << std::endl;
  }

  
  size_t count = 0;

  // store the line
  char str[1024];

  std::vector<std::vector<float>> vec(table.size());

  // read the remaining lines
  while (std::getline(rfile, line)) {
    
    strcpy(str, line.c_str());
    char *value_str = strtok(str, ",");
    std::vector<float> values(table.size());
    str[0] = '\0';
 
    size_t n = 0;
    while (value_str) {
      //table[std::next(table.begin(), n++)->first].push_back(atoi(value_str));
      vec[n++].push_back(atof(value_str));
      value_str = strtok(nullptr, ",");
    }

    if (count % 100000 == 0) {
      std::cerr << count << std::endl;
      }
    count++;
  }

  
  // move the values to the unordered_map
  for (auto& k : table) {
    auto it = std::find(elements.begin(), elements.end(), k.first);
    assert(it != elements.end());
    size_t index = std::distance(elements.begin(), it);
    k.second = vec[index];

    if (!markers.count(k.first) && k.first != x && k.first != y)
      meta.insert(k.first);
  }
  
  /* //debug print
  for (auto& k : markers)
    std::cerr << " M " << k << std::endl;
  for (auto& k : meta)
    std::cerr << " T " << k << std::endl;
  */
}

const std::vector<float>& CellTable::GetCoord(const char xyz) const {

  // check if this is expected
  if (xyz != 'x' && xyz != 'y' && xyz != 'z') {
    throw std::runtime_error("Coordinate is not equal to x, y, or z");
  }

  // Find the element with key "A" in the unordered_map
  std::string dim;
  switch (xyz) {
  case 'x' : dim = x; break;
  case 'y' : dim = y; break;
    //case 'z' : dim = z; break;    
  }
  auto it = table.find(dim);
  
  // Check if the element was found
  if (it == table.end()) {
    throw std::runtime_error("Element with key \"A\" not found in the unordered_map");
  }
  
  // Return a const reference to the float vector
  return it->second;
}



bool JsonReader::ReadData() {
  std::ifstream json_file(filename_);
  if (!json_file.is_open()) {
    std::cerr << "Error opening file " << filename_ << std::endl;
    return false;
  }
  
  Json::Value root;
  json_file >> root;
  
  x_ = root["x"].asString();
  y_ = root["y"].asString();
  z_ = root["z"].asString();    
  
  Json::Value markers = root["markers"];
  for (unsigned int i = 0; i < markers.size(); ++i) {
    markers_.insert(markers[i].asString());
  }
  
  Json::Value cellmeta = root["meta_cell"];
  for (unsigned int i = 0; i < cellmeta.size(); ++i) {
    cellmeta_.insert(cellmeta[i].asString());
  }
  
  return true;
}
