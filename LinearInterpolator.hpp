#include <map>
#include <vector>
 
 /**
  * Provides a basic interpolation mechanism in C++ using the STL.
  * Maybe not the fastest or most elegant method, but it works (for
  * linear interpolation!!), and is fast enough for a great deal of
  * purposes. It's also super easy to use, so that's a bonus.
  */
class LinearInterpolator {
	public:
		LinearInterpolator() {}
 
 		/**
 		 * Adds a data point to the interpolator for future interpolation
 		 * @param x The anchor point where this data is located
 		 * @param d A vector containing multiple "y" values to be interpolated (columns)
 		 */
		void addDataPoint(double x, std::vector<double> &d) {
			// just add it to the map
			data[x] = d;
		}
 
 		/**
 		 * Interpolates our data sequence at a given point, for a given column of data
 		 * @param  x      The anchor point to interpolate at
 		 * @param  column The column we want to interpolate on
 		 * @return y      The interpolated value for this column
 		 */
		double interpolate(double x, unsigned int column) {
			// loop through all the keys in the map
			// to find one that is greater than our intended value
			std::map<double, std::vector<double> >::iterator it = data.begin();
			bool found = false;
			while(it != data.end() && !found) {
				if(it->first >= x) {
					found = true;
					break;
				}
 
				// advance the iterator
				it++;
			}
 
			// check to see if we're outside the data range
			if(it == data.begin()) {
				return data.begin()->second[column];
			}
			else if(it == data.end()) {
				// move the point back one, as end() points past the list
				it--;
				return it->second[column];
			}
			// check to see if we landed on a given point
			else if(it->first == x) {
				return it->second[column];
			}
 
			// nope, we're in the range somewhere
			// collect some values
			double xb = it->first;
			double yb = it->second[column];
			it--;
			double xa = it->first;
			double ya = it->second[column];
 
			// and calculate the result!
			// formula from Wikipedia
			return (ya + (yb - ya) * (x - xa) / (xb - xa));
		}
 
	private:
		std::map<double, std::vector<double> > data;
};
