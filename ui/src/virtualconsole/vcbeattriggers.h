#pragma once
/*
  Q Light Controller Plus
  vcbeattriggers.h

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

#include "vcbutton.h"
#include "BTrack.h"

/** @addtogroup ui_vc_widgets
 * @{
 */

#define KXMLQLCVCBeatTriggers "BeatTriggers"

class VCBeatTriggers : public VCButton
{
    Q_OBJECT
    Q_DISABLE_COPY(VCBeatTriggers)

public:
    /*********************************************************************
     * Initialization
     *********************************************************************/
public:
    VCBeatTriggers(QWidget* parent, Doc* doc);
    ~VCBeatTriggers();

    /*********************************************************************
     * ID
     *********************************************************************/
public:
    /** @reimpl */
    void setID(quint32 id);

    /*********************************************************************
     * Clipboard
     *********************************************************************/
public:
    /** Create a copy of this widget to the given parent */
    VCWidget* createCopy(VCWidget* parent);

protected:
    /** Copy the contents for this widget from another widget */
    bool copyFrom(const VCWidget* widget);

    /*********************************************************************
     * audio capture
     *********************************************************************/
public:
    void enableCapture(bool enable);

    /** Method to toggle the enable button at a UI level.
     *  In this way we let Qt to handle the toggle signal and
     *  start the audio capture in the correct thread */
    //void toggleEnableButton(bool toggle);

public slots:
    //void slotEnableButtonToggled(bool toggle);

signals:
    void captureEnabled(bool enabled);

protected slots:
    void slotDetectBeat(const std::vector<double>& timeFrameData_p);

protected:
    AudioCapture*       m_inputCapture;
    BTrack              m_BTrack;
    std::vector<double> m_TimeFrameBuffer;
    std::vector<double> m_PrevTimeFrameBuffer;
    unsigned int        m_HopSize;
    unsigned int        m_FrameSize;

    /*********************************************************************
     * Properties
     *********************************************************************/
public:
    /** Edit this widget's properties */
    void editProperties();

    /*********************************************************************
     * Load & Save
     *********************************************************************/
public:
    /**
     * Load a VCBeatTriggers properties from an XML document node
     *
     * @param doc An XML document to load from
     * @param btn_root A VCBeatTriggers XML root node containing button properties
     * @return true if successful; otherwise false
     */
    bool loadXML(QXmlStreamReader &root);

    /**
     * Save a VCButton's properties to an XML document node
     *
     * @param doc The master XML document to save to
     * @param frame_root The button's VCFrame XML parent node to save to
     */
    bool saveXML(QXmlStreamWriter *doc);

};

/** @} */
