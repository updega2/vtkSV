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

#include "vtkSVRenderer.h"

#include "vtkCallbackCommand.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkTextActor.h"

#include "vtkSVIOUtils.h"
#include "vtkSVGeneralUtils.h"
#include "vtkSVGlobals.h"

#include <algorithm>

vtkStandardNewMacro(vtkSVRenderer);

vtkSVRenderer::vtkSVRenderer()
{
  this->WindowSize[0] = 800; this->WindowSize[1] = 600;
  this->WindowPosition[0] = 50; this->WindowPosition[1] = 50;
  this->Background[0] = 0.1; this->Background[1] = 0.1; this->Background[2] = 0.2;
  this->Annotations = 1;
  this->PointSmoothing = 1;
  this->LineSmoothing = 1;
  this->PolygonSmoothing = 0;
  this->TextInputMode = 0;
  this->ExitAfterTextInputMode = 1;
  this->ExitTextInputCallback = 0;
  this->InputPosition[0] = 0.25; this->InputPosition[1] = 0.1;
  this->Position[0] = 0.001; this->Position[1] = 0.05;

  this->Renderer = vtkRenderer::New();
  this->Renderer->SetBackground(this->Background);

  this->RenderWindow = vtkRenderWindow::New();
  this->RenderWindow->AddRenderer(this->Renderer);
  this->RenderWindow->SetSize(this->WindowSize[0],this->WindowSize[1]);
  this->RenderWindow->SetPosition(this->WindowPosition[0],this->WindowPosition[1]);
  this->RenderWindow->SetPointSmoothing(this->PointSmoothing);
  this->RenderWindow->SetLineSmoothing(this->LineSmoothing);
  this->RenderWindow->SetPolygonSmoothing(this->PolygonSmoothing);

  this->RenderWindowInteractor = vtkRenderWindowInteractor::New();

  this->RenderWindow->SetInteractor(this->RenderWindowInteractor);

  this->TrackballCamera = vtkInteractorStyleTrackballCamera::New();
  this->RenderWindowInteractor->SetInteractorStyle(this->TrackballCamera);
  this->RenderWindowInteractor->GetInteractorStyle()->KeyPressActivationOff();
  vtkNew(vtkCallbackCommand, charCallback);
  charCallback->SetCallback(vtkSVRenderer::CharCallback);
  this->RenderWindowInteractor->AddObserver("CharEvent", charCallback);
  vtkNew(vtkCallbackCommand, keyPressCallback);
  keyPressCallback->SetCallback(vtkSVRenderer::KeyPressCallback);
  this->RenderWindowInteractor->AddObserver("KeyPressEvent",keyPressCallback);

  vtkNew(vtkCallbackCommand, resetCameraCallback);
  resetCameraCallback->SetCallback(vtkSVRenderer::ResetCameraCallback);
  this->AddKeyBinding("r","Reset camera.",resetCameraCallback,"0");
  vtkNew(vtkCallbackCommand, quitRendererCallback);
  quitRendererCallback->SetCallback(vtkSVRenderer::QuitRendererCallback);
  this->AddKeyBinding("q","Quit renderer/proceed.",quitRendererCallback,"0");

  this->TextActor = vtkTextActor::New();
  this->TextActor->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
  this->TextActor->GetPosition2Coordinate()->SetCoordinateSystemToNormalizedViewport();
  this->TextActor->SetPosition(this->Position);
  this->Renderer->AddActor(this->TextActor);

  this->TextInputActor = vtkTextActor::New();
  this->TextInputActor->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
  this->TextInputActor->GetPosition2Coordinate()->SetCoordinateSystemToNormalizedViewport();
  this->TextInputActor->SetPosition(this->InputPosition);

}

vtkSVRenderer::~vtkSVRenderer()
{
  if (this->Renderer != NULL)
  {
    this->Renderer->Delete();
    this->Renderer = NULL;
  }
  if (this->RenderWindow != NULL)
  {
    this->RenderWindow->Delete();
    this->RenderWindow = NULL;
  }
  if (this->Renderer != NULL)
  {
    this->RenderWindowInteractor->Delete();
    this->RenderWindowInteractor = NULL;
  }
  if (this->TextActor != NULL)
  {
    this->TextActor->Delete();
    this->TextActor = NULL;
  }
  if (this->TextInputActor != NULL)
  {
    this->TextInputActor->Delete();
    this->TextInputActor = NULL;
  }
  if (this->TrackballCamera != NULL)
  {
    this->TrackballCamera->Delete();
    this->TrackballCamera = NULL;
  }
}

int vtkSVRenderer::Render(int interactive)
{
  if (interactive)
  {
    this->RenderWindowInteractor->Initialize();
  }
  this->RenderWindow->SetWindowName("vtkSV");

  // Somehow add all bindings
  //this->TextActor->SetInput('\n\n'.join(textActorInputsList));
  this->Renderer->AddActor(this->TextActor);

  this->RenderWindow->Render();

  if (interactive)
  {
    this->RenderWindowInteractor->Start();
  }

  return SV_OK;
}

int vtkSVRenderer::AddKeyBinding(std::string key, std::string text,
                                 vtkCallbackCommand *callback, std::string group)
{
  Binding newKey;
  newKey.key = key;
  newKey.text = text;
  newKey.callback = callback;
  newKey.group = group;

  this->KeyBindings.push_back(newKey);

  return SV_OK;
}

int vtkSVRenderer::RemoveKeyBinding(std::string key)
{
  std::vector<int> deleteList;
  for (int i=0; i<this->KeyBindings.size(); i++)
  {
    if (this->KeyBindings[i].key == key)
    {
      deleteList.push_back(i);
    }
  }

  for (int i=0; i<deleteList.size(); i++)
  {
      std::vector<Binding>::iterator delElem = this->KeyBindings.begin() + deleteList[i];
      this->KeyBindings.erase(delElem, this->KeyBindings.end());
  }

  return SV_OK;
}

int vtkSVRenderer::PromptAsync(std::string queryText, std::string callback)
{
  this->TextInputQuery = queryText;
  //this->ExitTextInputCallback = callback;
  this->UpdateTextInput();
  this->EnterTextInputMode(0);

  return SV_OK;
}

int vtkSVRenderer::EnterTextInputMode(int interactive)
{
  this->CurrentTextInput = "";
  this->Renderer->AddActor(this->TextInputActor);
  this->Renderer->RemoveActor(this->TextActor);
  this->UpdateTextInput();
  this->TextInputMode = 1;
  this->Render(interactive);

  return SV_OK;
}

int vtkSVRenderer::ExitTextInputMode()
{
  this->Renderer->RemoveActor(this->TextInputActor);
  this->Renderer->AddActor(this->TextActor);
  this->RenderWindow->Render();
  this->TextInputMode = 0;

  //if (this->ExitTextInputCallback)
  //  this->ExitTextInputCallback(this->CurrentTextInput);

  if (this->ExitAfterTextInputMode)
    this->RenderWindowInteractor->ExitCallback();

  return SV_OK;
}

int vtkSVRenderer::Close()
{
  //this->RenderWindowInteractor->Close();
}

void vtkSVRenderer::ResetCameraCallback( vtkObject* caller, long unsigned int vtkNotUsed(eventId), void* vtkNotUsed(clientData), void* vtkNotUsed(callData) )
{
  //this->Renderer->ResetCamera();
  //this->RenderWindow->Render();
}

void vtkSVRenderer::QuitRendererCallback( vtkObject* caller, long unsigned int vtkNotUsed(eventId), void* vtkNotUsed(clientData), void* vtkNotUsed(callData) )
{
  //this->PrintLog('Quit renderer')
  //this->Renderer->RemoveActor(this->TextActor)
  //this->RenderWindowInteractor->ExitCallback()
}

void vtkSVRenderer::KeyPressCallback( vtkObject* caller, long unsigned int vtkNotUsed(eventId), void* vtkNotUsed(clientData), void* vtkNotUsed(callData) )
{
  //std::cout << "Keypress callback" << std::endl;

  //vtkRenderWindowInteractor *iren =
  //  static_cast<vtkRenderWindowInteractor*>(caller);

  //std::cout << "Pressed: " << iren->GetKeySym() << std::endl;
  //std::string key = iren->GetKeySym();

  //if (key == "Escape")
  //{
  //  if (this->TextInputMode)
  //    this->TextInputMode = 0;
  //  else
  //    this->TextInputMode = 1;
  //}
  //if (this->TextInputMode)
  //{
  //  if (key == "Return" || key == "Enter")
  //  {
  //    this->ExitTextInputMode();
  //    return;
  //  }

  //  if (!strncmp(key, "KP_", 3))
  //    key = key.substr(3);

  //  if (key == "space")
  //    key = " ";
  //  else if (key == "minus" || key == "Subtract")
  //    key = "-";
  //  else if (key == "period" || key == "Decimal")
  //    key = ".";
  //  else if (key.length() > 1 && (key != "Backspace" || key != "BackSpace"))
  //    key = "";

  //  if (key != "Backspace" || key != "BackSpace")
  //  {
  //    std::string textInput = this->CurrentTextInput;
  //    if (textInput.length() > 0)
  //      this->CurrentTextInput = textInput;
  //  }
  //  else if (key != "")
  //    this->CurrentTextInput.append(key);

  //  this->UpdateTextInput();
  //  return;
  //}

  //int isKey = -1;
  //for (int i=0; i<this->KeyBindings.size(); i++)
  //{
  //  if (this->KeyBindings[i].key == key)
  //    isKey = i;
  //}

  //if (isKey != -1 && this->KeyBindings[isKey].callback != NULL)
  //{
  //  //this->KeyBindings[key]['callback'](obj);
  //  // SEt to something
  //}
  //else
  //{
  //  if (key == "plus")
  //    key = "+";
  //  if (key == "minus")
  //    key = "-";
  //  if (key == "equal")
  //    key = "=";
  //  if (isKey != -1 && this->KeyBindings[isKey].callback != NULL)
  //  {
  //    //this->KeyBindings[key]['callback'](obj);
  //    // SEt to something
  //  }
  //}
}

void vtkSVRenderer::UpdateTextInput()
{
  if (this->TextInputQuery != "")
  {
    std::string inputText = this->CurrentTextInput+"_";
    this->TextInputActor->SetInput(inputText.c_str());
    //this->TextInputActor->SetInput(this->TextInputQuery);
    this->Renderer->AddActor(this->TextInputActor);
  }
  else
  {
    this->Renderer->RemoveActor(this->TextInputActor);
  }
  this->RenderWindow->Render();
}

void vtkSVRenderer::CharCallback( vtkObject* caller, long unsigned int vtkNotUsed(eventId), void* vtkNotUsed(clientData), void* vtkNotUsed(callData) )
{
  return;
}

void vtkSVRenderer::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}
