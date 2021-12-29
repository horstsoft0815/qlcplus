/*
  Q Light Controller Plus
  vcbeattriggers.cpp

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

#include <QXmlStreamReader>
#include <QXmlStreamWriter>
#include <QDebug>

#include "vcbeattriggers.h"
#include "qlcfile.h"
#include "vcbeattriggersproperties.h"
#include "audiocapture.h"

/*****************************************************************************
 * Initialization
 *****************************************************************************/

VCBeatTriggers::VCBeatTriggers(QWidget* parent, Doc* doc)
    : VCButton(parent, doc)
    , m_inputCapture(NULL)
    , m_BTrack()
    , m_TimeFrameBuffer()
    , m_PrevTimeFrameBuffer()
    , m_HopSize(0)
    , m_FrameSize(0)
{
    setAction(Flash);
    setObjectName(VCBeatTriggers::staticMetaObject.className());
}

VCBeatTriggers::~VCBeatTriggers()
{
}

void VCBeatTriggers::setID(quint32 id)
{
    VCWidget::setID(id);

    if (caption().isEmpty())
        setCaption(tr("Beat Triggers %1").arg(id));
}

/*****************************************************************************
 * Clipboard
 *****************************************************************************/

VCWidget* VCBeatTriggers::createCopy(VCWidget* parent)
{
    Q_ASSERT(parent != NULL);

    VCBeatTriggers* button = new VCBeatTriggers(parent, m_doc);
    if (button->copyFrom(this) == false)
    {
        delete button;
        button = NULL;
    }

    return button;
}

bool VCBeatTriggers::copyFrom(const VCWidget* widget)
{
    const VCBeatTriggers* button = qobject_cast <const VCBeatTriggers*> (widget);
    if (button == NULL)
        return false;

    /* Copy button-specific stuff */
    setIconPath(button->iconPath());
    setKeySequence(button->keySequence());
    setFunction(button->function());
    enableStartupIntensity(button->isStartupIntensityEnabled());
    setStartupIntensity(button->startupIntensity());
    setAction(button->action());
    m_state = button->m_state;

    /* Copy common stuff */
    return VCWidget::copyFrom(widget);
}

/*********************************************************************
 * audio capture
 *********************************************************************/
void VCBeatTriggers::enableCapture(bool enable)
{
    // in case the audio input device has been changed in the meantime...
    QSharedPointer<AudioCapture> capture(m_doc->audioInputCapture());
    bool captureIsNew = m_inputCapture != capture.data();
    m_inputCapture = capture.data();

    if(!m_inputCapture)
    {
        return; // TODO: throw here
    }

    m_inputCapture->SetIsFrequencyAnalysisActive(false);
    m_inputCapture->SetIsTimeFrameAnalysisActive(true);

    m_FrameSize = m_inputCapture->GetBufferSize();
    m_HopSize = m_FrameSize/2;

    int hopSize = static_cast<int>(m_HopSize);
    int frameSize = static_cast<int>(m_FrameSize);
    m_BTrack = BTrack(hopSize, frameSize);
    m_TimeFrameBuffer.clear();
    m_PrevTimeFrameBuffer.clear();

    if (enable == true)
    {
        connect(m_inputCapture, SIGNAL(prepareTimeFrameData(const std::vector<double>&)),
                this, SLOT(slotDetectBeat(const std::vector<double>&)));

        this->blockSignals(true);
        //this->setChecked(true);
        this->blockSignals(false);

        emit captureEnabled(true);

        // Invalid ID: Stop every other widget
        emit functionStarting(Function::invalidId());
    }
    else
    {
        if (!captureIsNew)
        {
            disconnect(m_inputCapture, SIGNAL(prepareTimeFrameData(const std::vector<double>&)),
                       this, SLOT(slotDetectBeat(const std::vector<double>&)));
        }

        this->blockSignals(true);
        //this->setChecked(false);
        this->blockSignals(false);

        emit captureEnabled(false);
    }
}

/*void VCBeatTriggers::toggleEnableButton(bool toggle)
{
    if (mode() == Doc::Design)
        return;

    if (m_button)
        m_button->setChecked(toggle);
}

void VCBeatTriggers::slotEnableButtonToggled(bool toggle)
{
    if (mode() == Doc::Design)
        return;

    enableCapture(toggle);
}*/

void VCBeatTriggers::slotDetectBeat(const std::vector<double>& timeFrameData_p)
{
    // TODO: use local buffer to obtain overlap between old and new data --> hopsize

    if(m_TimeFrameBuffer.empty())
    {
        m_TimeFrameBuffer = timeFrameData_p;
        m_BTrack.processAudioFrame(m_TimeFrameBuffer.data());

        if(m_BTrack.beatDueInCurrentFrame())
        {
            pressFunction();
        }
    }
    else
    {
        std::vector<double> combinedFrame(m_TimeFrameBuffer);
        std::copy(timeFrameData_p.begin(), timeFrameData_p.end(), std::back_inserter(combinedFrame));

        std::vector<double> frameWithOverlap(m_FrameSize);
        std::size_t startPosition = m_HopSize;
        double* pStart = combinedFrame.data() + startPosition;
        double* pEnd = pStart + m_FrameSize;
        while(pEnd < (combinedFrame.data() + combinedFrame.size()))
        {
            std::copy(pStart, pEnd, frameWithOverlap.begin());

            m_BTrack.processAudioFrame(frameWithOverlap.data());

            if(m_BTrack.beatDueInCurrentFrame())
            {
                pressFunction();
            }

            startPosition += m_HopSize;
            pStart = combinedFrame.data() + startPosition;
            pEnd = pStart + m_FrameSize;
        }
    }



}

/*****************************************************************************
 * Properties
 *****************************************************************************/

void VCBeatTriggers::editProperties()
{
    VCBeatTriggersProperties prop(this, m_doc);
    if (prop.exec() == QDialog::Accepted)
        m_doc->setModified();
}

/*****************************************************************************
 * Load & Save
 *****************************************************************************/

bool VCBeatTriggers::loadXML(QXmlStreamReader &root)
{
    if (root.name() != KXMLQLCVCBeatTriggers)
    {
        qWarning() << Q_FUNC_INFO << "Beat Triggers node not found";
        return false;
    }

    return VCButton::loadXMLImpl(root);
}

bool VCBeatTriggers::saveXML(QXmlStreamWriter *doc)
{
    Q_ASSERT(doc != NULL);

    /* VC beat triggers entry */
    doc->writeStartElement(KXMLQLCVCBeatTriggers);

    return VCButton::saveXMLImpl(doc);
}
