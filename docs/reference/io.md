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

## cin
Standard C++ input.
```freefem
cin >> var;
```

## cout
Standard C++ ouptut.
```freefem
cout << "Some text" << endl;
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

## ifstream
Open a file in read mode.
```freefem
ifstream file("file.txt");
```

## ofstream
Open a file in write mode.
```freefem
ofstream file("file.txt");
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



