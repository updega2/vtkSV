/*=========================================================================
 *
 * Copyright (c) 2014 The Regents of the University of California.
 * All Rights Reserved.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *=========================================================================*/

#ifndef vtkSVRenderer_h
#define vtkSVRenderer_h

#include "vtkDataObject.h"
#include "vtkRenderer.h"
#include "vtkCallbackCommand.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkTextActor.h"
//#include "vtkvmtkComputationalGeometryWin32Header.h"
#include "vtkvmtkWin32Header.h"
#include "vtkSVCommonModule.h" // For exports

#include <vector>

struct Binding
{
  std::string key;
  std::string text;
  vtkCallbackCommand *callback = NULL;
  std::string group;
};

class VTKSVCOMMON_EXPORT vtkSVRenderer : public vtkDataObject
{
public:
  static vtkSVRenderer* New();

  vtkTypeMacro(vtkSVRenderer,vtkDataObject);


  int Render(int interactive);

  int AddKeyBinding(std::string key, std::string text,
                    vtkCallbackCommand *callback, std::string group);

  int RemoveKeyBinding(std::string key);

  int PromptAsync(std::string queryText, std::string callback);

  int EnterTextInputMode(int interactive);

  int ExitTextInputMode();

  int Close();

  void UpdateTextInput();

  static void ResetCameraCallback( vtkObject* caller, long unsigned int vtkNotUsed(eventId), void* vtkNotUsed(clientData), void* vtkNotUsed(callData) );

  static void QuitRendererCallback( vtkObject* caller, long unsigned int vtkNotUsed(eventId), void* vtkNotUsed(clientData), void* vtkNotUsed(callData) );

  static void KeyPressCallback( vtkObject* caller, long unsigned int vtkNotUsed(eventId), void* vtkNotUsed(clientData), void* vtkNotUsed(callData) );

  static void CharCallback( vtkObject* caller, long unsigned int vtkNotUsed(eventId), void* vtkNotUsed(clientData), void* vtkNotUsed(callData) );

  void PrintSelf(ostream& os, vtkIndent indent);

protected:
  vtkSVRenderer();
  ~vtkSVRenderer();

  vtkRenderer *Renderer;
  vtkRenderWindow *RenderWindow;
  vtkRenderWindowInteractor *RenderWindowInteractor;
  vtkInteractorStyleTrackballCamera *TrackballCamera;

  vtkTextActor *TextActor;
  vtkTextActor *TextInputActor;

  vtkCamera *Camera;

  int WindowSize[2];
  int WindowPosition[2];
  int Annotations;
  int PointSmoothing;
  int LineSmoothing;
  int PolygonSmoothing;
  int TextInputMode;
  int ExitAfterTextInputMode;
  int ExitTextInputCallback;
  int Interactive;

  double Background[3];
  double InputPosition[2];
  double Position[2];

  std::string TextInputQuery;
  std::string CurrentTextInput;

  std::vector<Binding> KeyBindings;

private:
  vtkSVRenderer(const vtkSVRenderer&);  // Not implemented.
  void operator=(const vtkSVRenderer&);  // Not implemented.
};

#endif
