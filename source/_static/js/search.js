const searchInput = document.getElementById('searchInput')
const searchResults = document.getElementById('searchResults')
let rootPath = ''

function initSearch(path) {
   rootPath = path
}

function search(event) {
   searchClean()

   if (!event.target)
      return

   if (!event.target.value)
      return

   const text = event.target.value

   if (text === '')
      return
   else
      searchLunr(text)
}

function searchClean() {
   searchResults.innerHTML = ''
}

function searchLunr(text) {
   const idx = lunr.Index.load(LUNR_DATA[0])
   // const results = idx.search('*'+text+'*')
   const results = idx.search(text)

   const resultsi = []
   results.forEach(function(result) {
      const ref = result['ref']
      const index = Number(ref)+1
      const idxi = lunr.Index.load(LUNR_DATA[index])
      // resultsi[index] = idxi.search('*'+text+'*')
      resultsi[index] = idxi.search(text)
   })

   const resultsHTML = parseLunrResults(results, resultsi, text)
   if (resultsHTML.length === 0)
      searchResults.innerHTML = '<p>No results</p>'
   else {
      for (let result of resultsHTML) {
         const div = document.createElement('div')
         div.className = 'search-result'

         const title = document.createElement('a')
         title.className = 'search-result-main-title'
         title.innerHTML = '<i class="fas fa-poll"></i>' + result.mainTitle
         title.href = result.link
         div.appendChild(title)

         for (let i = 0; i < result.titles.length; i++) {
            const subdiv = document.createElement('div')
            subdiv.className = 'search-result-sub'

            const title = document.createElement('a')
            title.className = 'search-result-title'
            title.href = result.links[i]

            const titleDiv = document.createElement('div')

            const titleText = document.createElement('p')
            titleText.className = 'search-result-title-text'
            titleText.innerHTML = result.titles[i]

            const preview = document.createElement('p')
            preview.className = 'search-result-preview'
            preview.innerHTML = result.previews[i]

            title.appendChild(titleDiv)
            titleDiv.appendChild(titleText)
            titleDiv.appendChild(preview)

            subdiv.appendChild(title)
            div.appendChild(subdiv)
         }

         searchResults.appendChild(div)
      }
      searchResults.style.display = 'block'
   }
}

function parseLunrResults(results, resultsi, text) {
   const html = []

   for (let i = 0; i < results.length; i++) {
      const id = results[i]['ref']
      const item = PREVIEW_LOOKUP[0][id]
      const mainTitle = item['t']
      const link = rootPath + item['l']

      const titlei = []
      const previewi = []
      const linki = []
      for (let k = 0; k < resultsi.length; k++) {
         const results = resultsi[k]
         if (results)
            for (var j = 0; j < results.length; j++) {
               const id = results[j]['ref']
               const item = PREVIEW_LOOKUP[k][id]
               const title = item['t']
               let preview = item['p']

               let lpreview = preview.toLowerCase()
               let ltext = text.toLowerCase()
               let index = lpreview.indexOf(ltext)
               preview = preview.slice(Math.max(0, index-124), Math.min(index+124, preview.length))

               lpreview = preview.toLowerCase()
               index = lpreview.indexOf(ltext)
               preview = preview.slice(0, index) + '<b>' + preview.slice(index, index+text.length) + '</b>' + preview.slice(index+text.length)

               const link = rootPath + item['l']
               titlei.push(title)
               previewi.push(preview)
               linki.push(link)
            }
      }

      html.push(
         {
            mainTitle: mainTitle,
            link: link,
            titles: titlei,
            previews: previewi,
            links: linki
         }
      )
   }

   return html
}

searchInput.addEventListener('input', function(event) {
   search(event)
})











// "use strict";
//
// var LUNR_CONFIG = {
//     "limit": 10,  // Max number of results to retrieve per page
//     "resultsElementId": "searchResults",  // Element to contain results
//     "countElementId": "resultCount"  // Element showing number of results
// };
//
// function search(text) {
//     clean()
//
//     if (text === '' || text === null)
//         return
//
//     searchLunr(text)
// }
//
// function clean() {
//     const countElement = document.getElementById(LUNR_CONFIG['countElementId'])
//     const resultsElement = document.getElementById(LUNR_CONFIG['resultsElementId'])
//
//     countElement.innerHTML = ''
//     resultsElement.innerHTML = ''
// }
//
// // Get URL arguments
// function getParameterByName(name, url) {
//     if (!url) url = 'window.location.href'
//     name = name.replace(/[\[\]]/g, "\\$&")
//     var regex = new RegExp("[?&]" + name + "(=([^&#]*)|&|#|$)"),
//         results = regex.exec(url)
//     if (!results) return null
//     if (!results[2]) return ""
//     return decodeURIComponent(results[2].replace(/\+/g, " "))
// }
//
// // Parse search results into HTML
// function parseLunrResults(results, query) {
//     var html = []
//     while(query.indexOf('*') !== -1)
//         query = query.replace('*', '');
//     for (var i = 0; i < results.length; i++) {
//         var id = results[i]["ref"]
//         var item = PREVIEW_LOOKUP[id]
//         var title = item["t"]
//         var preview = item["p"]
//         const pIndex = preview.toLowerCase().indexOf(query.toLowerCase())
//         preview = preview.slice(Math.max(pIndex, 64)-64, Math.max(pIndex, 128)+64)
//         var link = '../' + item["l"]
//         var result = ('<p><span class="result-title"><a href="' + link + '">'
//                     + title + '</a></span><br><span class="result-preview">'
//                     + preview + '</span></p>')
//         html.push(result)
//     }
//     if (html.length) {
//         return html.join("")
//     }
//     else {
//         return "<p>Your search returned no results.</p>";
//     }
// }
//
// function escapeHtml(unsafe) {
//     return unsafe
//         .replace(/&/g, "&amp;")
//         .replace(/</g, "&lt;")
//         .replace(/>/g, "&gt;")
//         .replace(/"/g, "&quot;")
//         .replace(/'/g, "&#039;");
// }
//
// function showResultCount(query, total, domElementId) {
//     if (total == 0) {
//         return
//     }
//
//     var s = "";
//     if (total > 1) {
//         s = "s"
//     }
//     var found = "<p>Found " + total + " result" + s
//     if (query != "" && query != null) {
//         query = escapeHtml(query)
//         var forQuery = ' for <span class="result-query">' + query + '</span>'
//     }
//     else {
//         var forQuery = ""
//     }
//     var element = document.getElementById(domElementId)
//     element.innerHTML = found + forQuery + "</p>"
// }
//
// function searchLunr(query) {
//     var idx = lunr.Index.load(LUNR_DATA)
//     // Write results to page
//     var results = idx.search(query)
//     var resultHtml = parseLunrResults(results, query)
//     var elementId = LUNR_CONFIG["resultsElementId"]
//     document.getElementById(elementId).innerHTML = resultHtml
//
//     var count = results.length;
//     showResultCount(query, count, LUNR_CONFIG["countElementId"])
// }
