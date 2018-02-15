## append
Append data to an existing file.
```ffpp
ofstream file("file.txt", append);
```

## binary
Write a file in binary.
```ffpp
ofstream file("file.btxt", binary);
```

## cin
Standard C++ input.
```ffpp
cin >> var;
```

## cout
Standard C++ ouptut.
```ffpp
cout << "Some text" << endl;
```

## flush
Flush the buffer of the file.
```ffpp
file.flush
```

## getline
Get the current line.
```ffpp
string s;
getline(file, s);
```

## ifstream
Open a file in read mode.
```
ifstream file("file.txt");
```

## ofstream
Open a file in write mode.
```ffpp
ofstream file("file.txt");
```

## seekg
Set the file position.
```ffpp
file.seekg(Pos);
```

## tellg
Get the file position.
```ffpp
int Pos = file.tellg();
```



