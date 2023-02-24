#ifndef CELL_HEADER2_H
#define CELL_HEADER2_H

#include <string>
#include <unordered_map>
#include <set>
#include <vector>

#include <iostream>
#include <fstream>

#include <json/json.h>

typedef std::set<std::string> Markers;
typedef std::set<std::string> Metas;

class CellTable {
    
public:

  CellTable() {}

  CellTable(const char* file, const char* markers_file);
  
  std::unordered_map<std::string, std::vector<float>> table;

  size_t num_cells() const { return table.at(x).size(); }

  using const_iterator = std::unordered_map<std::string, std::vector<float>>::const_iterator;
  const_iterator begin() const { return table.begin(); }
  const_iterator end() const { return table.end(); }

  // return a reference to the vector of x, y or z coordinates
  const std::vector<float>& GetCoord(const char xyz) const;

  // fill a vector of marker names
  // calling code is responsible for making a char* of same length as markers
  // calling code also needs to free each element and then the container char array
  void FillMarkerNames(char** array) const;

  // fill a vector of meta names
  // calling code is responsible for making a char* of same length as markers
  // calling code also needs to free each element and then the container char array
  void FillMetaNames(char** array) const;

  
  // return a vector of meta names
  std::vector<std::string> GetMetaNames() const;

  Metas meta;
  Markers markers;
  std::string x = "X_centroid";
  std::string y = "Y_centroid";

  // methods
  std::vector<double> XY() const;
  

};


class JsonReader {
public:
  JsonReader(const std::string& filename) : filename_(filename) {}
 
  std::string GetX() const { return x_; }
  std::string GetY() const { return y_; }
  const std::set<std::string>& GetMarkers() const { return markers_; }
  const std::set<std::string>& GetMetaCells() const { return cellmeta_; }  

  bool ReadData();
  
private:
  std::string filename_;
  std::string x_;
  std::string y_;
  std::string z_;  
  std::set<std::string> markers_;
  std::set<std::string> cellmeta_;  
};



#endif




