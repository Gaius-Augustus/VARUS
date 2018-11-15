/*
 * systemFunctions.cpp
 *
 *  Created on: 19.12.2017
 *      Author: willy
 */

#include <iostream>
#include <algorithm>
#include <string>
#include <dirent.h>

#include "../headers/debug.h"

using namespace std;

bool has_suffix(const string& s, const string& suffix)
{
    return (s.size() >= suffix.size()) && equal(suffix.rbegin(), suffix.rend(), s.rbegin());
}

vector<string> filesWithExtInFolder(const string& path, const string &ext){
	vector<string> files;
    DIR *dir = opendir(path.c_str());
    if(!dir)
    {
        DEBUG(0,"No files found with extension " << ext << " in folder " << path );
        return files;
    }

    dirent *entry;
    while((entry = readdir(dir))!=nullptr)
    {
        if(has_suffix(entry->d_name, ext))
        {
            string file = path + "/" + entry->d_name;
            files.push_back(file);
        }
    }

    closedir(dir);

    sort(files.begin(), files.end());
    return files;
}
