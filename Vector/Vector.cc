// Implementation of the templated Vector class
// ECE4893/8893 lab 3
// Amandine GOUT

#include <iostream> // debugging
#include "Vector.h"

// Your implementation here
// Fill in all the necessary functions below
using namespace std;

// Default constructor
template <typename T>
Vector<T>::Vector()
{
	count = 0;
	elements = NULL;
	reserved = 0;
}

// Copy constructor
template <typename T>
Vector<T>::Vector(const Vector& rhs)
{
	count = rhs.count;
	reserved = rhs.reserved;
	elements = (T*)malloc((rhs.reserved)*sizeof(T));
	for (size_t i=0; i<count;i++){
		new(&elements[i])T(rhs.elements[i]);
	}
}

// Assignment operator
template <typename T>
Vector<T>& Vector<T>::operator=(const Vector& rhs)
{
	// 3 cases : assigment to vector itself, vector not empty
	if (this != &rhs){
		if (elements != NULL){
			// call to destructor to delete current elements
			Clear();
			free(elements);
		}
		// Copy elements to new allocated memory
		count = rhs.count;
		reserved = rhs.reserved;
		elements = (T*)malloc((reserved)*sizeof(T));
		for (size_t i = 0; i<count;i++){
			new(&elements[i])T(rhs.elements[i]);
		}
	}
	return *this;
}

#ifdef GRAD_STUDENT
// Other constructors
template <typename T>
Vector<T>::Vector(size_t nReserved)
{ // Initialize with reserved memory
	count = 0;
	reserved = nReserved;
	elements = (T*)malloc(nReserved*sizeof(T));
}

template <typename T>
Vector<T>::Vector(size_t n, const T& t)
{ // Initialize with "n" copies of "t"
	elements = (T*)malloc(n*sizeof(T));
	reserved = n;
	for (size_t i=0;i<n;i++){
		new(&elements[i])T(t);
	}
	count = n;
}

template <typename T>
void Vector<T>::Reserve(size_t n)
{ // Reserve extra memory
	if (elements==NULL){
		elements = (T*)malloc((reserved+n)*sizeof(T));
	} else {
		T* newelements = (T*)malloc((reserved+n)*sizeof(T));
		if (count>0){
			for (size_t i=0;i<count;i++){
				new(&newelements[i])T(elements[i]);
				elements[i].~T();
			}
		}
		free(elements);
		elements = newelements;
	}
	reserved = reserved+n;
}

#endif

// Destructor
template <typename T>
Vector<T>::~Vector()
{
	Clear();
	free(elements);
}

// Add and access front and back
template <typename T>
void Vector<T>::Push_Back(const T& rhs)
{
	// case when enough memory
	if (reserved>count){
		new(&elements[count])T(rhs);
	} else { // case when not enough memory
		T* newelements = (T*)malloc((count+1)*sizeof(T));
		reserved = count+1;
		for (size_t i=0; i<count; i++){
			new(&newelements[i])T(elements[i]);
			elements[i].~T();
		}
		new(&newelements[count])T(rhs);
		free(elements);
		elements = newelements;
	}
	count = count+1;
}

template <typename T>
void Vector<T>::Push_Front(const T& rhs)
{
	// case when enough memory
	if (reserved>count){
		for (size_t i=1; i<count+1; i++){
			new(&elements[count-i+1])T(elements[count-i]);
			T(elements[count-i]).~T();
		}
		new(&elements[0])T(rhs);
		count = count+1;
	} else { // case when not enough memory
		T* newelements = (T*)malloc((count+1)*sizeof(T));
		reserved = count+1;
		new(&newelements[0])T(rhs);
		for (size_t i=1; i<count+1; i++){
			new(&newelements[i])T(elements[i-1]);
			elements[i-1].~T();
		}
		free(elements);
		elements = newelements;
	}
	count = count+1;
}

template <typename T>
void Vector<T>::Pop_Back()
{ // Remove last element
	if (count == 0){
		cout << "Empty vector: cannot pop back element " << endl;
		abort();
	} else {
		elements[count-1].~T();
		count = count -1;
	}
}

template <typename T>
void Vector<T>::Pop_Front()
{ // Remove first element
	if (count == 0){
		cout << "Empty vector: cannot pop front element " << endl;
		abort();
	} else {
		for (size_t i=0; i<count-1; i++){
			elements[i].~T();
			new(&elements[i])T(elements[i+1]);
		}
		elements[count-1].~T();
		count = count -1;
	}
}

// Element Access
template <typename T>
T& Vector<T>::Front() const
{
	if (elements == NULL){
		cout << "Front access: NULL pointer!" << endl;
		abort();
	} else {
		return (*(elements));
	}
}

// Element Access
template <typename T>
T& Vector<T>::Back() const
{
	if (elements == NULL){
		cout << "Back access : NULL pointer!" << endl;
		abort();
	} else {
		return (*(elements+count-1));
	}
	
}

template <typename T>
const T& Vector<T>::operator[](size_t i) const
{ // const element access
	if (i < 0 || i >= count)
    {
    cout << "ERROR INDEXING!" << endl;
    abort();
    }
    return (*(elements+i));
}

template <typename T>
T& Vector<T>::operator[](size_t i)
{//nonconst element access
	if (i < 0 || i >= count)
    {
    cout << "ERROR INDEXING!" << endl;
    abort();
    }
    return (*(elements+i));
}

template <typename T>
size_t Vector<T>::Size() const
{
	return count;
}

template <typename T>
bool Vector<T>::Empty() const
{
	if (count == 0){
		return true;
	}
	else{ 
		return false;
	}
		
}

// Implement clear
template <typename T>
void Vector<T>::Clear()
{
	reserved = count;
	for (size_t i=0;i<count;++i){
		elements[count-i-1].~T();
	}
	count = 0;
}

// Iterator access functions
template <typename T>
VectorIterator<T> Vector<T>::Begin() const
{
  return VectorIterator<T>(elements);
}

template <typename T>
VectorIterator<T> Vector<T>::End() const
{
  return VectorIterator<T>(elements+count);
}

#ifdef GRAD_STUDENT
// Erase and insert
template <typename T>
void Vector<T>::Erase(const VectorIterator<T>& it)
{
	if (elements != NULL) {
		count = count -1;
		VectorIterator<T> it1 = this->Begin();
		bool half = false;
		size_t i = 0;
		VectorIterator<T> it2 = this->End();
		while (it1 != it2){
			if (it1 == it){
				half = true;
			} 
			if (half == true) {
				elements[i].~T();
				new(&elements[i])T(elements[i+1]);
			}
			it1++;
			i++;
		}
		elements[count].~T();
	} else {
		cout << "NULL pointer ! Cannot Erase element" << endl;
		abort();
	}
}

template <typename T>
void Vector<T>::Insert(const T& rhs, const VectorIterator<T>& it)
{
	if (elements != NULL) {
		size_t i = 0;
		bool half = false;
		VectorIterator<T> it1 = this->Begin();
		T* newelements = (T*)malloc((count+1)*sizeof(T));
		reserved  = count+1;
		while (it1 != this->End()){
			if (it1 == it){
				new(&newelements[i])T(rhs);
				half = true;
			} 
			if (half == true) {
				new(&newelements[i+1])T(elements[i]);
				elements[i].~T();
			} else {
				new(&newelements[i])T(elements[i]);
				elements[i].~T();
			}
			it1++;
			i++;
		}
		free(elements);
		elements = newelements;
		count = count + 1;
	} else {
		cout << "NULL pointer ! Cannot Erase element" << endl;
		abort();
	}
}
#endif

// Implement the iterators

// Constructors
template <typename T>
VectorIterator<T>::VectorIterator()
{
	current = NULL;
}

template <typename T>
VectorIterator<T>::VectorIterator(T* c)
{
	current = c;
}

// Copy constructor
template <typename T>
VectorIterator<T>::VectorIterator(const VectorIterator<T>& rhs)
{
	current = rhs.current;
}

// Iterator dereferencing operator
template <typename T>
T& VectorIterator<T>::operator*() const
{
	return *current;
}

// Prefix increment : call ++it
template <typename T>
VectorIterator<T>  VectorIterator<T>::operator++()
{
	current = current+1;
	return *this;
}

// Postfix increment : call it++
template <typename T>
VectorIterator<T> VectorIterator<T>::operator++(int)
{
	// Postfix returns the original value prior to increment
	VectorIterator<T> previous = VectorIterator(current);
	//VectorIterator<T> previous = *this;
	current = current+1;
	return previous;
}

// Comparison operators
template <typename T>
bool VectorIterator<T>::operator !=(const VectorIterator<T>& rhs) const
{
	if (current != rhs.current){
		return true;
	} else {
		return false;
	}
}

template <typename T>
bool VectorIterator<T>::operator ==(const VectorIterator<T>& rhs) const
{
	if (current == rhs.current){
		return true;
	} else {
		return false;
	}
}




