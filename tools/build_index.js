const path = require('path')
const fs = require('fs')
const lunr = require('lunr')
const cheerio = require('cheerio')

const HTML_FOLDER = '../build/html'
const SEARCH_FIELDS = ['title', 'body']
const EXCLUDE_FILES = ['search.html', 'genindex.html']
const OUTPUT_INDEX = 'lunr_index'

function isHtml(filename) {
   lower = filename.toLowerCase()
   return (lower.endsWith('.htm') || lower.endsWith('.html'))
}

function findHtml(folder) {
   if (!fs.existsSync(folder)) {
      console.log('Could not find folder: ', folder)
      return
   }

   const files = fs.readdirSync(folder)
   const htmls = []
   for (let i = 0; i < files.length; i++) {
      const filename = path.join(folder, files[i])
      const stat = fs.lstatSync(filename)
      if (stat.isDirectory()) {
         const recursed = findHtml(filename)
         for (var j = 0; j < recursed.length; j++) {
            recursed[j] = path.join(files[i], recursed[j]).replace(/\\/g, '/')
         }
         htmls.push.apply(htmls, recursed)
      }
      else if (isHtml(filename) && !EXCLUDE_FILES.includes(files[i])) {
         htmls.push(files[i])
      }
   }

   return htmls
}

function readHtml(root, file, fileId) {
   const filename = path.join(root, file)
   const txt = fs.readFileSync(filename).toString()

   const $ = cheerio.load(txt)

   function removeTocTree(element) {
      element.children().filter(function(i, el) {
         const className = $(this).attr('class')
         if (className && className.includes('toctree-wrapper'))
            $(this).remove()
         else
            removeTocTree($(this))
      })
   }
   removeTocTree($('#content'))

   let title = $('title').text()
   if (typeof title == 'undefined') title = file

   let body = $('#content').text()
   if (typeof body == 'undefined') body = ''

   const data = [{
      'id': fileId,
      'link': file,
      't': title,
      'b': body
   }]

   return data
}

function parseAnchor(text) {
   text = text.toLowerCase()
   text = text.replace(' ', '-')
   text = text.replace(/[&\/\\#,+()$~%.'":*?<>{}]/g,'-')
   return text
}

function readSingleHtml(root, file) {
   const filename = path.join(root, file)
   const txt = fs.readFileSync(filename).toString()

   const $ = cheerio.load(txt)

   function removeTocTree(element) {
      element.children().filter(function(i, el) {
         const className = $(this).attr('class')
         if (className && className.includes('toctree-wrapper'))
            $(this).remove()
         else
            removeTocTree($(this))
      })
   }
   removeTocTree($('#content'))

   const data = []

   const content = $('#content')
   const sections = content.find('.section')

   sections.each(function (i, elem) {
      $(this).find('.section').remove()
   })

   sections.each(function (i, elem) {
      let title = $(this).find('h1')
      if (title.length === 0)
         title = $(this).find('h2')
         if (title.length === 0)
            title = $(this).find('h3')
            if (title.length === 0)
               title = $(this).find('h4')
               if (title.length === 0)
                  title = $(this).find('h5')
                  if (title.length === 0)
                     title = $(this).find('h6')

      title = title.text()
      const body = $(this).text()

      data.push({
         'id': i,
         'link': file + '#' + parseAnchor(title),
         't': title,
         'b': body
      })
   })

   return data
}

function buildIndex(docs) {
   const idx = lunr(function () {
      this.ref('id')
      for (let i = 0; i < SEARCH_FIELDS.length; i++) {
         this.field(SEARCH_FIELDS[i].slice(0, 1))
      }
      docs.forEach(function (doc) {
         this.add(doc)
      }, this)
   })
   return idx
}

function buildPreviews(docs) {
   const result = {}
   for (let i = 0; i < docs.length; i++) {
      const doc = docs[i]
      const preview = doc['b']
      result[doc['id']] = {
         't': doc['t'],
         'p': preview,
         'l': doc['link']
      }
   }
   return result
}

function main() {
   files = findHtml(HTML_FOLDER)
   let docs = []
   const idxi = []
   const previewsi = []
   console.log('Building index for these files:')
   for (let i = 0; i < files.length; i++) {
      console.log('    ' + files[i])
      docs = docs.concat(readHtml(HTML_FOLDER, files[i], i))

      const doc = readSingleHtml(HTML_FOLDER, files[i])
      idxi.push(buildIndex(doc))
      previewsi.push(buildPreviews(doc))
   }
   const idx = buildIndex(docs)
   const previews = buildPreviews(docs)
   let js = 'const LUNR_DATA = [' + JSON.stringify(idx)
   for (let i = 0; i < idxi.length; i++)
      js += ',\n' + JSON.stringify(idxi[i])
   js += '];\n'
   js += 'const PREVIEW_LOOKUP = [' + JSON.stringify(previews)
   for (let i = 0; i < previewsi.length; i++)
      js += ',\n' + JSON.stringify(previewsi[i])
   js += '];'
   fs.writeFile(OUTPUT_INDEX+'.js', js, function(err) {
      if(err) {
         return console.log(err)
      }
      console.log('Index saved as ' + OUTPUT_INDEX)
   })
}

main()
