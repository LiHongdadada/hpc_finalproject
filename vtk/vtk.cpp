#include <hdf5.h>
#include "vtkIntArray.h"
#include "vtkDoubleArray.h"
#include "vtkPolyData.h"
#include "vtkCellData.h"
#include "vtkCellArray.h"
#include "vtkPointData.h"
#include "vtkTetra.h"
#include "vtkGenericCell.h"
#include "vtkUnstructuredGrid.h"
#include "vtkUnstructuredGridWriter.h"
#include "vtkXMLUnstructuredGridWriter.h"
#include "vtkXMLUnstructuredGridReader.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkXMLGenericDataObjectReader.h"
#include <vector>
#include <stdio.h>
using namespace std;
#define FILE "explicit.h5"
void add_int_PointData(vtkPointSet *const &grid_w,
                       const std::vector<int> &ptdata, const std::string &dataname);

void add_int_CellData(vtkPointSet *const &grid_w,
                      const std::vector<int> &cldata, const std::string &dataname);
void add_double_PointData(vtkPointSet *const &grid_w,
                         const std::vector<double> &ptdata, const std::string &dataname);
int main(int argc, char **argv)
{
    hid_t file_id, dataset_T_id, dataset_times_id; /* identifiers */
    double times[3];
    file_id = H5Fopen(FILE, H5F_ACC_RDWR, H5P_DEFAULT);
    dataset_times_id = H5Dopen(file_id, "/iteration_times", H5P_DEFAULT);
    H5Dread(dataset_times_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, times);
    float h = times[1];
    int n = 1 / h;
    int num_of_nodes = (n + 1) * (n + 1), num_of_elements = n * n;
    double T[num_of_nodes];
    float nodes[num_of_nodes][3];
    int elements[num_of_elements][4];

    for (int i = 0; i < num_of_nodes; i++)
    {
        int a = i / (n + 1);
        nodes[i][0] = (i - a * (n + 1)) * h;
        nodes[i][1] = a * h;
        nodes[i][2] = 0;
    }
    for (int i = 0; i < num_of_elements; i++)
    {
        int a = i / n;
        elements[i][0] = i + a;
        elements[i][1] = i + a + 1;
        elements[i][2] = i + a + n + 2;
        elements[i][3] = i + a + n + 1;
    }
    dataset_T_id = H5Dopen(file_id, "/T", H5P_DEFAULT);
    H5Dread(dataset_T_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, T);

    const std::string filename("temperature");
    vtkUnstructuredGrid *grid_w = vtkUnstructuredGrid::New();
    vtkPoints *ppt = vtkPoints::New();
    ppt->SetDataTypeToDouble();

    for (int i = 0; i < num_of_nodes; i++)
    {
        ppt->InsertPoint(i, nodes[i][0], nodes[i][1], nodes[i][2]);
    }
    grid_w->SetPoints(ppt);
    ppt->Delete();

    vtkCell *cl = vtkTetra::New();
    for (int ii = 0; ii < num_of_elements; ++ii)
    {

        for (int k = 0; k < 4; k++)
        {
            cl->GetPointIds()->SetId(k, elements[ii][k]);
        }
        grid_w->InsertNextCell(cl->GetCellType(), cl->GetPointIds());
    }
    cl->Delete();

    vector<int> node_id;
    for (int i = 0; i < num_of_nodes; i++)
    {
        node_id.push_back(i);
    }

    add_int_PointData(grid_w, node_id, "GlobalNodeID");

    vector<int> cell_id;
    for (int i = 0; i < num_of_elements; i++)
    {
        cell_id.push_back(i);
    }
    add_int_CellData(grid_w, cell_id, "GlobalCellID");
    vector<double> Tem(num_of_nodes);
    memcpy(&Tem[0], T, sizeof(T));
    add_double_PointData(grid_w, Tem, "Temperature");

    vtkUnstructuredGridWriter *writer = vtkUnstructuredGridWriter::New();
    std::string name_to_write(filename);

    name_to_write.append(".vtk");
    writer->SetFileName(name_to_write.c_str());

    writer->SetInputData(grid_w);
    writer->Write();
    writer->Delete();

    grid_w->Delete();
    H5Dclose(dataset_T_id);
    H5Dclose(dataset_times_id);
    H5Fclose(file_id);
    return 0;
}
void add_int_PointData(vtkPointSet *const &grid_w,
                       const std::vector<int> &ptdata, const std::string &dataname)
{
    vtkIntArray *data = vtkIntArray::New();
    data->SetNumberOfComponents(1);
    data->SetName(dataname.c_str());

    for (unsigned int ii = 0; ii < ptdata.size(); ++ii)
        data->InsertComponent(ii, 0, ptdata[ii]);

    grid_w->GetPointData()->AddArray(data);
    data->Delete();
}
void add_double_PointData(vtkPointSet *const &grid_w,
                         const std::vector<double> &ptdata, const std::string &dataname)
{
    vtkDoubleArray *data = vtkDoubleArray::New();
    data->SetNumberOfComponents(1);
    data->SetName(dataname.c_str());

    for (unsigned int ii = 0; ii < ptdata.size(); ++ii)
        data->InsertComponent(ii, 0, ptdata[ii]);

    grid_w->GetPointData()->AddArray(data);
    data->Delete();
}
void add_int_CellData(vtkPointSet *const &grid_w,
                      const std::vector<int> &cldata, const std::string &dataname)
{
    vtkIntArray *data = vtkIntArray::New();
    data->SetNumberOfComponents(1);
    data->SetName(dataname.c_str());

    for (unsigned int ii = 0; ii < cldata.size(); ++ii)
        data->InsertComponent(ii, 0, cldata[ii]);

    grid_w->GetCellData()->AddArray(data);
    data->Delete();
}
