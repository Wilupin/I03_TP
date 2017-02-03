#include "plot.hpp"
#include <cstdio>

#ifdef NO_GRAPHICS
Plot::Plot()
{
}

Plot::~Plot()
{
}

void Plot::save()
{
}

void Plot::init(size_t , size_t , size_t , Poisson_Parameters & )
{
}

void Plot::terminate() {
}

void Plot::process(const Values &)
{
}

#else

#include <vtkVersion.h>
#include <vtkSmartPointer.h>

#include <vtkImageData.h>
#include <vtkImageActor.h>
#include <vtkCamera.h>
#include <vtkImageMapper.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkMath.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkImageMapToColors.h>
#include <vtkImageMapper3D.h>
#include <vtkPolyDataMapper.h>
#include <vtkImageMapToColors.h>
#include <vtkImageProperty.h>
#include <vtkInteractorStyleImage.h>
#include <vtkLookupTable.h>
#include <vtkSmartPointer.h>
#include <vtkVersion.h>

struct PlotInternal {
  vtkSmartPointer<vtkImageData> grid;
  vtkSmartPointer<vtkRenderer> renderer; 
  vtkSmartPointer<vtkRenderWindow> renderWindow;
  vtkSmartPointer<vtkLookupTable> lookupTable;
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor;
  double Mmin, Mmax;  
};

Plot::Plot() {
  m_plot = new PlotInternal;
  m_plot->Mmin = 0.0;
  m_plot->Mmax = 0.0;
  m_order = 0;
}

Plot::~Plot() {
  delete m_plot;
}

void Plot::init(size_t n, size_t m, size_t p, Poisson_Parameters & P)
{
  m_resultPath = P.resultPath();

  vtkSmartPointer<vtkImageData> grid =
    vtkSmartPointer<vtkImageData>::New();
  m_plot->grid = grid;

  int imageExtent[6] = { 0, int(n)-1, 0, int(m)-1, 0, int(p/2) - 1 };
  grid->SetExtent(imageExtent);

#if VTK_MAJOR_VERSION <= 5
  grid->SetNumberOfScalarComponents(1);
  grid->SetScalarTypeToDouble();
#else
  grid->AllocateScalars(VTK_DOUBLE,1);
#endif

  size_t i,j,k;
  for (k = 0; k < p/2; k++)
    for (j = 0; j < m; j++)
    for (i = 0; i < n; i++) {
      double* pixel = static_cast<double*>(grid->GetScalarPointer(i,j,k));
      pixel[0] = 0.0;
    }
  
  // Map the scalar values in the image to colors with a lookup table:
  vtkSmartPointer<vtkLookupTable> lookupTable
    = vtkSmartPointer<vtkLookupTable>::New();
  m_plot->lookupTable = lookupTable;
  
  lookupTable->SetNumberOfTableValues(64);
  lookupTable->SetRange(0.0, 1.0);
  lookupTable->SetTableValue( 0, 0.00000,   0.00000,   0.50000, 0.8);
  lookupTable->SetTableValue( 1, 0.00000,   0.00000,   0.56349, 0.8);
  lookupTable->SetTableValue( 2, 0.00000,   0.00000,   0.62698, 0.8);
  lookupTable->SetTableValue( 3, 0.00000,   0.00000,   0.69048, 0.8);
  lookupTable->SetTableValue( 4, 0.00000,   0.00000,   0.75397, 0.8);
  lookupTable->SetTableValue( 5, 0.00000,   0.00000,   0.81746, 0.8);
  lookupTable->SetTableValue( 6, 0.00000,   0.00000,   0.88095, 0.8);
  lookupTable->SetTableValue( 7, 0.00000,   0.00000,   0.94444, 0.8);
  lookupTable->SetTableValue( 8, 0.00000,   0.00794,   1.00000, 0.8);
  lookupTable->SetTableValue( 9, 0.00000,   0.07143,   1.00000, 0.8);
  lookupTable->SetTableValue(10, 0.00000,   0.13492,   1.00000, 0.8);
  lookupTable->SetTableValue(11, 0.00000,   0.19841,   1.00000, 0.8);
  lookupTable->SetTableValue(12, 0.00000,   0.26190,   1.00000, 0.8);
  lookupTable->SetTableValue(13, 0.00000,   0.32540,   1.00000, 0.8);
  lookupTable->SetTableValue(14, 0.00000,   0.38889,   1.00000, 0.8);
  lookupTable->SetTableValue(15, 0.00000,   0.45238,   1.00000, 0.8);
  lookupTable->SetTableValue(16, 0.00000,   0.51587,   1.00000, 0.8);
  lookupTable->SetTableValue(17, 0.00000,   0.57937,   1.00000, 0.8);
  lookupTable->SetTableValue(18, 0.00000,   0.64286,   1.00000, 0.8);
  lookupTable->SetTableValue(19, 0.00000,   0.70635,   1.00000, 0.8);
  lookupTable->SetTableValue(20, 0.00000,   0.76984,   1.00000, 0.8);
  lookupTable->SetTableValue(21, 0.00000,   0.83333,   1.00000, 0.8);
  lookupTable->SetTableValue(22, 0.00000,   0.89683,   1.00000, 0.8);
  lookupTable->SetTableValue(23, 0.00000,   0.96032,   1.00000, 0.8);
  lookupTable->SetTableValue(24, 0.02381,   1.00000,   0.97619, 0.8);
  lookupTable->SetTableValue(25, 0.08730,   1.00000,   0.91270, 0.8);
  lookupTable->SetTableValue(26, 0.15079,   1.00000,   0.84921, 0.8);
  lookupTable->SetTableValue(27, 0.21429,   1.00000,   0.78571, 0.8);
  lookupTable->SetTableValue(28, 0.27778,   1.00000,   0.72222, 0.8);
  lookupTable->SetTableValue(29, 0.34127,   1.00000,   0.65873, 0.8);
  lookupTable->SetTableValue(30, 0.40476,   1.00000,   0.59524, 0.8);
  lookupTable->SetTableValue(31, 0.46825,   1.00000,   0.53175, 0.8);
  lookupTable->SetTableValue(32, 0.53175,   1.00000,   0.46825, 0.8);
  lookupTable->SetTableValue(33, 0.59524,   1.00000,   0.40476, 0.8);
  lookupTable->SetTableValue(34, 0.65873,   1.00000,   0.34127, 0.8);
  lookupTable->SetTableValue(35, 0.72222,   1.00000,   0.27778, 0.8);
  lookupTable->SetTableValue(36, 0.78571,   1.00000,   0.21429, 0.8);
  lookupTable->SetTableValue(37, 0.84921,   1.00000,   0.15079, 0.8);
  lookupTable->SetTableValue(38, 0.91270,   1.00000,   0.08730, 0.8);
  lookupTable->SetTableValue(39, 0.97619,   1.00000,   0.02381, 0.8);
  lookupTable->SetTableValue(40, 1.00000,   0.96032,   0.00000, 0.8);
  lookupTable->SetTableValue(41, 1.00000,   0.89683,   0.00000, 0.8);
  lookupTable->SetTableValue(42, 1.00000,   0.83333,   0.00000, 0.8);
  lookupTable->SetTableValue(43, 1.00000,   0.76984,   0.00000, 0.8);
  lookupTable->SetTableValue(44, 1.00000,   0.70635,   0.00000, 0.8);
  lookupTable->SetTableValue(45, 1.00000,   0.64286,   0.00000, 0.8);
  lookupTable->SetTableValue(46, 1.00000,   0.57937,   0.00000, 0.8);
  lookupTable->SetTableValue(47, 1.00000,   0.51587,   0.00000, 0.8);
  lookupTable->SetTableValue(48, 1.00000,   0.45238,   0.00000, 0.8);
  lookupTable->SetTableValue(49, 1.00000,   0.38889,   0.00000, 0.8);
  lookupTable->SetTableValue(50, 1.00000,   0.32540,   0.00000, 0.8);
  lookupTable->SetTableValue(51, 1.00000,   0.26190,   0.00000, 0.8);
  lookupTable->SetTableValue(52, 1.00000,   0.19841,   0.00000, 0.8);
  lookupTable->SetTableValue(53, 1.00000,   0.13492,   0.00000, 0.8);
  lookupTable->SetTableValue(54, 1.00000,   0.07143,   0.00000, 0.8);
  lookupTable->SetTableValue(55, 1.00000,   0.00794,   0.00000, 0.8);
  lookupTable->SetTableValue(56, 0.94444,   0.00000,   0.00000, 0.8);
  lookupTable->SetTableValue(57, 0.88095,   0.00000,   0.00000, 0.8);
  lookupTable->SetTableValue(58, 0.81746,   0.00000,   0.00000, 0.8);
  lookupTable->SetTableValue(59, 0.75397,   0.00000,   0.00000, 0.8);
  lookupTable->SetTableValue(60, 0.69048,   0.00000,   0.00000, 0.8);
  lookupTable->SetTableValue(61, 0.62698,   0.00000,   0.00000, 0.8);
  lookupTable->SetTableValue(62, 0.56349,   0.00000,   0.00000, 0.8);
  lookupTable->SetTableValue(63, 0.50000,   0.00000,   0.00000, 0.8);

  lookupTable->Build();

  // Pass the original image and the lookup table to a filter to create
  // a color image:
  vtkSmartPointer<vtkImageMapToColors> scalarValuesToColors =
    vtkSmartPointer<vtkImageMapToColors>::New();
  scalarValuesToColors->SetLookupTable(lookupTable);
  scalarValuesToColors->PassAlphaToOutputOn();
#if VTK_MAJOR_VERSION <= 5
  scalarValuesToColors->SetInput(grid);
#else
  scalarValuesToColors->SetInputData(grid);
#endif
   
  // Create an image actor
  vtkSmartPointer<vtkImageActor> imageActor =
    vtkSmartPointer<vtkImageActor>::New();
  imageActor->GetMapper()->SetInputConnection
    (scalarValuesToColors->GetOutputPort());
  imageActor->GetProperty()->SetInterpolationTypeToNearest();

  // Visualize
  vtkSmartPointer<vtkRenderer> renderer =
    vtkSmartPointer<vtkRenderer>::New();
  m_plot->renderer = renderer;

  renderer->AddActor(imageActor);
  
  vtkSmartPointer<vtkRenderWindow> renderWindow =
    vtkSmartPointer<vtkRenderWindow>::New();
  m_plot->renderWindow = renderWindow;
  
  renderWindow->SetSize(600,600);
  renderWindow->AddRenderer(renderer); 
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  m_plot->renderWindowInteractor = renderWindowInteractor;

  vtkSmartPointer<vtkInteractorStyleImage> style =
    vtkSmartPointer<vtkInteractorStyleImage>::New();
  
  renderWindowInteractor->SetInteractorStyle(style);
  renderWindowInteractor->SetRenderWindow(renderWindow);
  renderWindowInteractor->Initialize();
}

void Plot::terminate() {
  m_plot->renderWindowInteractor->Start();
}

void Plot::process(const Values &u)
{  
  size_t i,j,k, n = u.size(0), m = u.size(1), p = u.size(2); 
  double Mmin_new = 0.0, Mmax_new = 0.0;
  
  for (i = 0; i < n; i++)
    for (j = 0; j < m; j++)
      for (k = 0; k < p/2; k++) {
	double v = u(i,j,k+p/2);
	if (v < Mmin_new) Mmin_new = v;
	else if (v > Mmax_new) Mmax_new = v;
	
	double* pixel
	  = static_cast<double*>(m_plot->grid->GetScalarPointer(i,j,k));
	pixel[0] = v;
    }

  if (Mmax_new > m_plot->Mmax * 2 || Mmax_new < m_plot->Mmax / 2)
    m_plot->Mmax = Mmax_new;
  if (Mmin_new < m_plot->Mmin * 2 || Mmin_new < m_plot->Mmin / 2)
    m_plot->Mmin = Mmin_new;
  
  m_plot->lookupTable->SetRange(m_plot->Mmin, m_plot->Mmax);
  m_plot->lookupTable->Build();


  // Pass the original image and the lookup table to a filter to create
  // a color image:
  vtkSmartPointer<vtkImageMapToColors> scalarValuesToColors =
    vtkSmartPointer<vtkImageMapToColors>::New();
  scalarValuesToColors->SetLookupTable(m_plot->lookupTable);
  scalarValuesToColors->PassAlphaToOutputOn();
#if VTK_MAJOR_VERSION <= 5
  scalarValuesToColors->SetInput(m_plot->grid);
#else
  scalarValuesToColors->SetInputData(m_plot->grid);
#endif

  m_plot->grid->Modified();
  m_plot->renderWindow->Render();

  save();
  m_order++;
}

#include <sstream>
void Plot::save() {

  std::ostringstream s;
  s << m_resultPath << "/plot_"
    << m_order << ".vti";

  vtkSmartPointer<vtkXMLImageDataWriter> writer =
    vtkSmartPointer<vtkXMLImageDataWriter>::New();
  writer->SetFileName(s.str().c_str());
  writer->SetDataModeToAscii();

#if VTK_MAJOR_VERSION <= 5
  writer->SetInputConnection(m_plot->grid->GetProducerPort());
#else
  writer->SetInputData(m_plot->grid);
#endif
  writer->Write();
}

#endif
