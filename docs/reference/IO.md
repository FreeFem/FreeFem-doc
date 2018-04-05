See [I/O example](../example/#i/o)

See [File stream example](../examples/file-stream).

## cout
Standard C++ output device (default: console).

```freefem
cout << "Some text" << endl;
```

## cin
Standard C++ input device (default: keyboard).

```freefem
cin >> var;
```

## endl

End of line.

```freefem
cout << "Some text" << endl;
```

## ifstream
Open a file in read mode.
```freefem
ifstream file("file.txt");
```

!!!note
	A file is closed at the end of a block.

## ofstream
Open a file in write mode.
```freefem
ofstream file("file.txt");
```

!!!note
	A file is closed at the end of a block.

## append
Append data to an existing file.
```freefem
ofstream file("file.txt", append);
```

## binary
Write a file in binary.
```freefem
ofstream file("file.btxt", binary);
```

## seekg
Set the file position.
```freefem
file.seekg(Pos);
```

## tellg
Get the file position.
```freefem
int Pos = file.tellg();
```

## flush
Flush the buffer of the file.
```freefem
file.flush
```

## getline
Get the current line.
```freefem
string s;
getline(file, s);
```

## Output format

In the descriptions below, `f` is an output stream, for example `cout` or a `ofstream`.

All this methods, excepted the first, return a stream, so they can be chained:
```freefem
cout.scientific.showpos << 3 << endl;
```

### precision

Set the number of digits printed to the right of the decimal point. This applies to all subsequent floating point numbers written to that output stream. However, this won't make floating-point "integers" print with a decimal point. It's necessary to use `:::freefem fixed` for that effect.

```freefem
int np = f.precision(n)
```

### scientific

Formats floating-point numbers in scientific notation

```freefem
f.scientific
```

### fixed

Used fixed point notation for floating-point numbers. Opposite of scientific.

```freefem
f.fixed
```

### showbase

Converts insertions to an external form that can be read according to the `C++` lexical conventions for integral constants. By default, showbase is not set.

```freefem
f.showbase
```

### noshowbase

Unset `:::freefem showbase` flags.

```freefem
f.noshowbase
```

### showpos

Inserts a plus sign (+) into a decimal conversion of a positive integral value.

```freefem
f.showpos
```

### noshowpos

Unset `:::freefem showpos` flags.

```freefem
f.noshowpos
```

### default

Reset all the previous flags to the default expect precision.

```freefem
f.default
```

### setw

Behaves as if member width were called with `n` as argument on the stream on which it is inserted as a manipulator (it can be inserted on output streams).

```freefem
f.setw(n)
```
