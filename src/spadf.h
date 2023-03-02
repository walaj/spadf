#ifndef SPADF_H
#define SPADF_H

#include "cell2.h"
#include <netcdf.h>

///////
// MACROS
///////
#define NC_DEFINE_DIM(name, len, id)			\
  if (len > 0) {							\
    if (nc_def_dim(ncid, name, len, &id) != NC_NOERR) {			\
      std::cerr << "Error defining " << name << " dimension: " << nc_strerror(ncerr) << std::endl; \
      return;								\
    }									\
  }

#define NC_DEFINE_VAR(name, type, rank, dimids, id) \
if (nc_def_var(ncid, name, type, rank, dimids, &id) != NC_NOERR) \
{ \
  std::cerr << "Error defining " << name << " variable: " << nc_strerror(ncerr) << std::endl; \
  return; \
}

#define NC_WRITE_VAR_FLOAT(name, id, data) \
if (nc_put_var_float(ncid, id, data) != NC_NOERR) \
{ \
  std::cerr << "Error writing " << name << " data: " << nc_strerror(ncerr) << std::endl; \
}

#define NC_WRITE_VAR_STRING(name, id, data) \
if (nc_put_var_string(ncid, id, data) != NC_NOERR) \
{ \
  std::cerr << "Error writing " << name << " data: " << nc_strerror(ncerr) << std::endl; \
}


class spadf {

 public:

  spadf(){}

  // a null-terminated string containing file name and extension
  spadf(const char* file, const CellTable& table);

  // add the 2D marker fluorescence table
  int AddMarkerTable(const CellTable& table);

  int AddMetaTable(const CellTable& table);
  
  // add the 2D cell meta table
  int AddCellMetaTable(const CellTable& table);

  // close the file
  int Close();

  // generate the knn graph in marker space
  int knn_marker();

  // subsample the cells
  int subsample_cells(int n, int seed);

  // return rows of marker data at a time
  // maker sure that the data is allocated by the calling function
  int get_marker_rows(size_t row, size_t nrows, float* data);
  
 private:

  int ncid = -1;
  
  // data dimids
  int markers_dimid = -1;
  int cells_dimid = -1;
  int coords_dimid = -1;

  // data lens
  size_t markers_len = 0;
  size_t cells_len = 0;
  size_t coords_len = 0;
  
  // meta dimids
  int meta_markers_dimid = -1;
  int meta_cells_dimid = -1;

  // meta lens
  size_t meta_markers_len = 0;
  size_t meta_cells_len = 0;
  
  // 1D: data varids (meta data)
  int markers_varid = -1;
  int cells_varid = -1;
  int meta_cells_varid = -1;  
  
  // 2D: data varids (fluorescence and spatial)
  int marker_table_varid = -1;
  int coord_table_varid = -1;
  
  // 2D: cells meta varids
  int cells_meta_table_varid = -1;

  //////
  // methods
  //////

  // checks that the input table is the expected shape and size
  // and returns the number of cells
  int _table_check(const CellTable& table);
 
  
};

#endif
