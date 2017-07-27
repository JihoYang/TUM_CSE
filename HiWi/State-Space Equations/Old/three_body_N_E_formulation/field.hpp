#ifndef FIELD_H_
#define FIELD_H_

//system includes
#include <cstdlib>
#include <cassert>

template<typename T>
class field;

template<typename T>
class field {

private:
	size_t _width;
	size_t _height;
	size_t _size;
	T* _data;
public:
	field(){
	}
	//constructor with specific size.
	field(size_t w, size_t h) :
			_width(w), _height(h), _size(w * h) {
		_data = new T[_size];
	}

	//~field() { delete _data; }

	inline const size_t getWidth() const {
		return _width;
	}

	inline const size_t getHeight() const {
		return _height;
	}

	inline const size_t getSize() const {
		return _size;
	}

	//access element (x,y) (2D Access)
	inline T & operator()(size_t x, size_t y) {

		assert(x<_width);
		assert(y<_height);

		return _data[y * _width + x];

	}

	//access element [x] (1D-Access)
	inline T & operator[](size_t x) {

		assert(x<_size);

		return _data[x];

	}

	//access element (x,y) (2D Access)
	inline const T operator()(size_t x, size_t y) const {
		//assert(x<_width);
		//assert(y<_height);
		if (x >= _width)
			x = _width - 1;
		if (y >= _height)
			y = _height - 1;

		return _data[y * _width + x];

	}

	//access element [x] (1D-Access)
	inline const T operator[](size_t x) const {

		assert(x<_size);

		return _data[x];

	}
};

#endif // FIELD_H_
