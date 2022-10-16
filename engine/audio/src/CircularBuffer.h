//=======================================================================
/** @file CircularBuffer.h
 *  @brief A class for calculating onset detection functions
 *  @author Adam Stark
 *  @copyright Copyright (C) 2008-2014  Queen Mary University of London
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
//=======================================================================

#ifndef CircularBuffer_h
#define CircularBuffer_h

#include <vector>

//=======================================================================
/** A circular buffer that allows you to add new samples to the end
 * whilst removing them from the beginning. This is implemented in an
 * efficient way which doesn't involve any memory allocation
 */
class CircularBuffer
{
public:
    
    /** Constructor */
    CircularBuffer()
     :  writeIndex (0)
    {
        
    }
    
    /** Access the ith element in the buffer */
    double& operator[] (int i)
    {
        int index = (i + writeIndex) % buffer.size();
        return buffer[index];
    }
    
    const double& operator[] (int i) const
    {
        int index = (i + writeIndex) % buffer.size();
        return buffer[index];
    }

    /** Add a new sample to the end of the buffer */
    void addSampleToEnd (double v)
    {
        buffer[writeIndex] = v;
        writeIndex = (writeIndex + 1) % buffer.size();

        if(currentSize < buffer.size())
        {
            ++currentSize;
        }
    }
    
    /** Resize the buffer */
    void resize (int size)
    {
        buffer.resize (size);
        writeIndex = 0;
        currentSize = size>0 ? 1 : 0;
    }

    std::size_t getSize() const { return buffer.size(); }
    std::size_t getCurrentSize() const { return currentSize; }
    int getWriteIndex() const { return writeIndex; }

    const std::vector<double>& getBuffer() const { return buffer; }
private:
    
    std::vector<double> buffer;
    std::size_t currentSize;
    int writeIndex;
};

namespace CircularBufferUtils
{
struct LinearFit
{
    double slope = 0.;
    double intercept = 0.;
};

std::vector<double> GetDifferenceVec(const CircularBuffer& buffer_p);
LinearFit GetLinearFit(const double timeStep_p, const CircularBuffer& buffer_p);
}


#endif /* CircularBuffer_hpp */
