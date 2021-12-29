/*
  Q Light Controller
  vcbuttonproperties.cpp

  Copyright (c) Sebastian Moeckel

  Licensed under the Apache License, Version 2.0 (the "License");
  you may not use this file except in compliance with the License.
  You may obtain a copy of the License at

      http://www.apache.org/licenses/LICENSE-2.0.txt

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.
*/


#include "vcbeattriggersproperties.h"
#include "inputselectionwidget.h"

VCBeatTriggersProperties::VCBeatTriggersProperties(
        VCBeatTriggers* button, Doc* doc)
    : VCButtonProperties(button, doc)
{
    // hide different button modes, beat triggers only supports "flash"
    m_toggle->hide();
    m_flash->setChecked(true);
    m_flash->hide();
    m_blackout->hide();
    m_stopAll->hide();
    groupBox_4->hide();

    // hide external input, beat triggers uses audio source as input
    m_inputSelWidget->hide();

    this->adjustSize();

    setWindowTitle(QApplication::translate("VCBeatTriggers", "Beat Triggers properties", nullptr));
}

VCBeatTriggersProperties::~VCBeatTriggersProperties()
{
}
